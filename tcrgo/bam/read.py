import pysam
import textwrap
import re
from collections import defaultdict, Counter, OrderedDict
from .cdr3 import Transcript, CDR3
from .reference import ReferenceDict
import tcrgo.bio as bio
from Bio.Seq import Seq
from typing import List, Dict, Iterator, Set, Tuple, Union, Optional, Any
from ..log import Log

log = Log("root")
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

class Read(object):
	def __init__(self, query_name: str):
		self.query_name = query_name
		self.alignments = list()
		self.top_J = None
		self.top_J_reference = None
		self.top_V = None
		self.top_V_reference = None
		self.cdr3 = None
		self.has_region = defaultdict(bool)
		self.ties = Counter()

	def __len__(self) -> int:
		return len(self.alignments)

	def __iter__(self) -> Iterator:
		return iter(self.alignments)

	def __getitem__(self, item: str):
		return getattr(self, item)

	def __setitem__(self, item: str, value: Any):
		setattr(self, item, value)

	def append(self, alignment: AlignedSegment):
		self.alignments.append(alignment)

	def compare_metrics(self, top_alignment: AlignedSegment, alignment: AlignedSegment) \
		-> Tuple[Optional[AlignedSegment], Optional[AlignedSegment], Optional[str]]:
		"""Attempt to break a score tie by comparing other alignment metrics"""
		winner = None
		loser = None
		method = None
		if alignment.get_tag("NM") != top_alignment.get_tag("NM"): # Edit Distance
			method = "EditDistance"
			winner, loser = sorted((top_alignment, alignment), key=lambda a: a.get_tag("NM"))
		elif alignment.get_tag("UQ") != top_alignment.get_tag("UQ"): # Phred Quality
			method = "Phred"
			loser, winner = sorted((top_alignment, alignment), key=lambda a: a.get_tag("UQ"))
		elif alignment.query_alignment_length != top_alignment.query_alignment_length:
			method = "AlignmentLength"
			loser, winner = sorted((top_alignment, alignment), key=lambda a: a.query_alignment_length)
		return winner, loser, method

	def break_tie(self, top_alignment: AlignedSegment, alignment: AlignedSegment) \
		-> Union[List[AlignedSegment], AlignedSegment]:
		"""Handle a alignment score tie between the top_alignment and the contending alignment"""
		winner, loser, method = self.compare_metrics(top_alignment, alignment)
		if winner is not None:
			self.ties[(winner.reference_name, loser.reference_name, method)] += 1
		else:
			winner = sorted([top_alignment, alignment], key=lambda a: a.reference_name)
		return winner

	def break_tie_multi(self, top_alignments: List[AlignedSegment], alignment: AlignedSegment) \
		-> Union[List[AlignedSegment], AlignedSegment]:
		"""Handle a alignment score tie when top_alignment is already an unresolved tie."""
		top_alignment_rep = top_alignments[0]
		winner, loser, method = self.compare_metrics(top_alignment_rep, alignment)
		if winner is not None:
			if winner.reference_name == top_alignment_rep.reference_name: # winner is top_alignments
				winner = top_alignments
				ref_winner = tuple([top_alignment.reference_name for top_alignment in top_alignments])
				ref_loser = loser.reference_name
			else: # winner is the singular alignment
				ref_winner = winner.reference_name
				ref_loser = tuple([top_alignment.reference_name for top_alignment in top_alignments])
			self.ties[(ref_winner, ref_loser, method)] += 1
		else:
			top_alignments.append(alignment)
			winner = sorted(top_alignments, key=lambda a: a.reference_name)
		return winner

	def compare_scores(self, top_alignment: Union[AlignedSegment, List[AlignedSegment]], \
		top_score: int, alignment: AlignedSegment) \
		-> Union[Tuple[AlignedSegment, int], Tuple[List[AlignedSegment], int]]:
		"""Compare alignment scores and return the best alignment(s) and the new best score."""
		score = alignment.get_tag("AS")
		if score >= top_score:
			if score == top_score:
				if isinstance(top_alignment, list):
					top_alignment = self.break_tie_multi(top_alignment, alignment)
				else:
					top_alignment = self.break_tie(top_alignment, alignment)
			else:
				top_alignment = alignment
				top_score = score
		return top_alignment, top_score

	def parse_alignments(self, bam_indexed: IndexedReads, refdict: ReferenceDict):
		"""Find the top V and J alignments, count segment identities"""
		top_score_J = 0
		top_score_V = 0
		for alignment in bam_indexed.find(self.query_name):
			self.alignments.append(alignment)
			#self.unique_subregions.add(alignment.reference_name)
			for v in ("TRAV", "TRBV"):
				if v in alignment.reference_name:
					self.has_region[v] = True
					if alignment.reference_end < len(refdict[alignment.reference_name].seq_nt) - 20:
						continue
					self.top_V, top_score_V = self.compare_scores(self.top_V, top_score_V, alignment)
					break
			else:
				for j in ("TRAJ", "TRBJ"):
					if j in alignment.reference_name:
						self.has_region[j] = True
						self.top_J, top_score_J = self.compare_scores(self.top_J, top_score_J, alignment)
						break
				else:
					for c in ("TRAC", "TRBC"):
						if c in alignment.reference_name:
							self.has_region[c] = True
							break
					else:
						self.has_region["UNKN"] = True

	def get_deletion_positions(self, cigar: str) -> Dict[int, int]:
		deletion_positions = dict()
		num = ""
		position = 0
		for char in cigar:
			if  47 < ord(char) < 58: # If number 0-9
				num += char
			elif char == 'D':
				deletion_positions[position] = int(num) # Key: start, value: width
			else: # Letter that isn't D
				position += int(num) - 1
				num = ""
		return deletion_positions

	def construct_transcript(self) -> Transcript:
		query_start = self.top_V.reference_start - self.top_V.query_alignment_start
		abstract_V = self.top_V_reference.seq_nt[:query_start]
		abstract_J = self.top_J_reference.seq_nt[self.top_J.reference_end:]
		transcript = Transcript(
			seq_nt = abstract_V + Seq(self.top_V.query_sequence) + abstract_J,
			query_start_nt = query_start,
			V_end_aa = (query_start + self.top_V.query_alignment_end) // 3,
			J_start_aa = (query_start + self.top_J.query_alignment_start) // 3,
			orf = self.top_V_reference.orf
		)
		transcript.get_translation()
		return transcript

	def get_cdr3(self):
		log.verbose(self.query_name)
		transcript = self.construct_transcript()
		query_cdr3_start_nt = self.top_V.query_alignment_start \
			- self.top_V.reference_start + self.top_V_reference.cdr3_positions
		cdr3_start_aa = (transcript.query_start_nt + query_cdr3_start_nt) // 3
		if transcript.seq_aa[cdr3_start_aa] != 'C':
			cdr3_starts_aa = transcript.get_cdr3_starts()
			if cdr3_starts_aa:
				cdr3_start_aa = cdr3_starts_aa[0] # min
		cdr3_start_nt = cdr3_start_aa * 3 + transcript.orf
		
		cdr3_ends_aa = list()
		if self.top_J_reference.cdr3_positions[transcript.orf]:
			query_cdr3_ends_nt = [
				self.top_J.query_alignment_start - self.top_J.reference_start + end \
					for end in self.top_J_reference.cdr3_positions[transcript.orf]
			]
			for end in query_cdr3_ends_nt:
				end_nt = transcript.query_start_nt + end
				end_aa = end_nt // 3
				if transcript.seq_aa[end_aa] == 'F':
					cdr3_ends_aa.append(end_aa)
		if not cdr3_ends_aa:
			cdr3_ends_aa = transcript.get_cdr3_ends()
			if not cdr3_ends_aa:
				cdr3_ends_aa = [cdr3_start_aa + 15]
		cdr3_end_nt = cdr3_ends_aa[-1] * 3 + transcript.orf

		cdr3 = CDR3(
			transcript.seq_nt[cdr3_start_nt:cdr3_end_nt+3],
			cdr3_start_nt, cdr3_end_nt+2, transcript
		)
		cdr3.get_translation()
		if not cdr3.is_in_frame():
			log.error("Out of frame!")
		if len(transcript.get_stops()) == 0:
			cdr3.is_productive = True
		self.cdr3 = cdr3
		log.verbose(self.VJ_alignment_stats(), indent=0)
		log.verbose(self.VJ_alignment_diagram_full(), indent=0)

	def VJ_alignment_diagram(self) -> str:
		"""Includes just the part of the reference that aligns to the query"""
		ref_seq_V = self.top_V_reference.seq_nt
		ref_seq_J = self.top_J_reference.seq_nt
		query_cdr3_start = self.cdr3.query_start
		query_cdr3_end = self.cdr3.query_end
		return textwrap.dedent(
			f"""
			{' ' * self.top_V.query_alignment_start}{ref_seq_V[self.top_V.reference_start:self.top_V.reference_end]}
			{self.visualize_cigar(self.top_V.cigarstring)}
			{' ' * self.top_V.query_alignment_start}{'?'*(self.top_V.reference_end-self.top_V.reference_start)}
			{' ' * query_cdr3_start}*
			{self.top_V.query_sequence}
			{' ' * query_cdr3_end}*
			{' ' * self.top_J.query_alignment_start}{'?'*(self.top_J.reference_end-self.top_J.reference_start)}
			{self.visualize_cigar(self.top_J.cigarstring)}
			{' ' * self.top_J.query_alignment_start}{ref_seq_J[self.top_J.reference_start:self.top_J.reference_end]}
			"""
		)

	def VJ_alignment_diagram_full(self) -> str:
		"""
		Contains the entire V and J references. It is recommended you look at the output
		in an editor that has text-wrapping disabled.
		"""
		query_start = self.top_V.reference_start - self.top_V.query_alignment_start
		j_start = query_start + self.top_J.query_alignment_start - self.top_J.reference_start
		transcript_translations = bio.get_frame_translations(self.cdr3.transcript.seq_nt)
		return textwrap.dedent(
			f"""
			{self.top_V_reference.seq_nt}\t{self.top_V_reference.name}
			{' ' * query_start}{bio.visualize_alignment(self.top_V)}\t{self.top_V.cigarstring}\t{self.top_V.get_tag("MD")}
			{' ' * self.cdr3.start}S
			{' ' * query_start}{self.top_V.query_sequence}\tQUERY {self.query_name}
			{' ' * self.cdr3.start}{str(self.cdr3.seq_nt)}\t\t\t{str(self.cdr3.seq_aa)}
			{' ' * self.cdr3.end}E
			{' ' * query_start}{bio.visualize_alignment(self.top_J)}\t{self.top_J.cigarstring}\t{self.top_J.get_tag("MD")}
			{' ' * j_start}{self.top_J_reference.seq_nt}\t{self.top_J_reference.name}
			{self.cdr3.transcript.seq_nt}\tCONSTRUCTED TRANSCRIPT V->J
			{bio.visualize_aminoacids(transcript_translations[0])}
			{' ' + bio.visualize_aminoacids(transcript_translations[1])}
			{"  " + bio.visualize_aminoacids(transcript_translations[2])}
			"""
		)
		# {' ' * self.top_V.reference_start}{self.top_V_reference.seq_nt[self.top_V.reference_start:self.top_V.reference_end]}
		# {' ' * query_start}{' '*self.top_J.query_alignment_start}{self.top_J_reference.seq_nt[self.top_J.reference_start:self.top_J.reference_end]}

	@log.fdoc
	def VJ_alignment_stats(self) -> str:
		def center(text: Optional[str]) -> str:
			if text is None:
				text = "N/A"
			return "{:^34}".format(text)
		return \
			f"""
			{'_'*80}
			| RNAME |{center(self.top_V.reference_name)}||{center(self.top_J.reference_name)}|
			| RBEG  |{center(self.top_V.reference_start)}||{center(self.top_J.reference_start)}|
			| REND  |{center(self.top_V.reference_end)}||{center(self.top_J.reference_end)}|
			| RLEN  |{center(str(len(self.top_V_reference.seq_nt)))}||{center(str(len(self.top_J_reference.seq_nt)))}|
			| RCDR3 |{center(str(self.top_V_reference.cdr3_positions))}||{center(str(self.top_J_reference.cdr3_positions))}|
			|{' '*78}|
			| QBEG  |{center(self.top_V.query_alignment_start)}||{center(self.top_J.query_alignment_start)}|
			| QEND  |{center(self.top_V.query_alignment_end)}||{center(self.top_J.query_alignment_end)}|
			| QALEN |{center(self.top_V.query_alignment_length)}||{center(self.top_J.query_alignment_length)}|
			|{' '*78}|
			| CIGAR |{center(self.top_V.cigarstring)}||{center(self.top_J.cigarstring)}|
			|  MD   |{center(self.top_V.get_tag('MD'))}||{center(self.top_J.get_tag('MD'))}|
			| SCORE |{center(self.top_V.get_tag('AS'))}||{center(self.top_J.get_tag('AS'))}|
			| EDIST |{center(self.top_V.get_tag('NM'))}||{center(self.top_J.get_tag('NM'))}|
			| PHRED |{center(self.top_V.get_tag('UQ'))}||{center(self.top_J.get_tag('UQ'))}|
			{'_'*80}
			"""
	
	@log.fdoc
	def alignment_info(self, alignment: AlignedSegment) -> str:
		return \
			f"""
			QNAME:	{alignment.query_name}
			REF:	{alignment.reference_name}
			SCORE:	{alignment.get_tag("AS")}
			CIGAR:	{alignment.cigarstring}
			PHRED:	{alignment.get_tag("UQ")}
			EDIST:	{alignment.get_tag("NM")}
			FLAG:	{alignment.flag}
			MAPQ:	{alignment.mapping_quality}
			"""

	def alignment_info_line(self, alignment: AlignedSegment) -> str:
		if alignment == None:
			return f"NONE{'	'*7}"
		return '\t'.join([
			alignment.query_name,
			alignment.reference_name,
			str(alignment.get_tag("AS")),
			alignment.cigarstring,
			str(alignment.get_tag("UQ")),
			str(alignment.get_tag("NM")),
			str(alignment.flag),
			str(alignment.mapping_quality)
		]) + '\n'

	@log.fdoc
	def __str__(self) -> str:
		return \
			f"""
			QNAME:	{self.query_name}
			NALIGN:	{len(self)}
			TOPV:	{self.top_V.reference_name if self.top_V is not None else self.top_V}
			TOPJ:	{self.top_J.reference_name if self.top_J is not None else self.top_J}
			"""
