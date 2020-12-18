import pysam
import textwrap
import re
from collections import defaultdict, Counter, OrderedDict
from .cdr3 import Transcript, CDR3
from .reference import ReferenceDict
import tcrgo.bio as bio
from Bio.Seq import Seq
# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set, Tuple, Union, Optional, Any
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

from ..log import Log
log = Log(name=__name__)
log.proceed()

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
		log.info(f"transcript.query_start_nt={transcript.query_start_nt} transcript.V_end_aa={transcript.V_end_aa} transcript.J_start_aa={transcript.J_start_aa} transcript.orf={transcript.orf}")
		return transcript

	def get_cdr3(self):
		log.info(self.query_name)
		log.info(self.VJ_alignment_stats())
		transcript = self.construct_transcript()
		query_cdr3_start_nt = self.top_V.query_alignment_start \
			- self.top_V.reference_start + self.top_V_reference.cdr3_positions
		log.verbose(f"query_cdr3_start_nt={query_cdr3_start_nt}")
		cdr3_start_aa = (transcript.query_start_nt + query_cdr3_start_nt) // 3
		log.verbose(f"cdr3_start_aa={cdr3_start_aa}")

		log.verbose(f"transcript.seq_nt={transcript.seq_nt} {len(transcript.seq_nt)}")
		log.verbose(f"transcript.seq_aa={transcript.seq_aa} {len(transcript.seq_aa)}")
		log.verbose(f"transcript.seq_aa[cdr3_start_aa]={transcript.seq_aa[cdr3_start_aa]}")
		if transcript.seq_aa[cdr3_start_aa] != 'C':
			cdr3_starts_aa = transcript.get_cdr3_starts()
			log.verbose(f"cdr3_start_aa={cdr3_starts_aa}")
			if cdr3_starts_aa:
				cdr3_start_aa = cdr3_starts_aa[0] # min
		cdr3_start_nt = cdr3_start_aa * 3 + transcript.orf
		log.verbose(f"cdr3_start_aa={cdr3_start_nt}")
		
		cdr3_ends_aa = list()
		log.verbose(f"self.top_J_reference.cdr3_positions[transcript.orf]={self.top_J_reference.cdr3_positions[transcript.orf]}")
		if self.top_J_reference.cdr3_positions[transcript.orf]:
			query_cdr3_ends_nt = [
				self.top_J.query_alignment_start - self.top_J.reference_start + end \
					for end in self.top_J_reference.cdr3_positions[transcript.orf]
			]
			log.verbose(f"query_cdr3_ends_nt={query_cdr3_ends_nt}")
			for end in query_cdr3_ends_nt:
				end_nt = transcript.query_start_nt + end
				end_aa = end_nt // 3
				if transcript.seq_aa[end_aa] == 'F':
					cdr3_ends_aa.append(end_aa)
		log.verbose(f"cdr3_ends_aa={cdr3_ends_aa}")
		if not cdr3_ends_aa:
			cdr3_ends_aa = transcript.get_cdr3_ends()
			log.verbose(f"cdr3_ends_aa={cdr3_ends_aa}")
			if not cdr3_ends_aa:
				cdr3_ends_aa = [cdr3_start_aa + 15]
		cdr3_end_nt = cdr3_ends_aa[-1] * 3 + transcript.orf
		log.verbose(f"cdr3_end_nt={cdr3_end_nt}")

		log.verbose(f"len(transcript.seq_nt)={len(transcript.seq_nt)}")
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
		print(self.VJ_alignment_diagram_full())

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
		def center(text: str) -> str:
			return "{:^34}".format(text)
		#if self.top_V_reference.cdr3_positions is not None:
		#	query_cdr3_starts = str([self.top_V.query_alignment_start - self.top_V.reference_start + start \
		#		for start in self.top_V_reference.cdr3_positions])
		#else:
		query_cdr3_starts = "N/A"
		#if self.top_J_reference.cdr3_positions is not None:
		#	query_cdr3_ends = str([self.top_J.query_alignment_start - self.top_J.reference_start + end \
		#		for end in self.top_J_reference.cdr3_positions])
		#else:
		query_cdr3_ends = "N/A"
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
			| QCDR3 |{center(query_cdr3_starts)}||{center(query_cdr3_ends)}|
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

#########################################################################
#	Below is unused or deprecated code
#########################################################################

def score_cdr3_candidates(self, query_cdr3_starts: List[int], query_cdr3_ends: List[int]) -> List[CDR3]:
	cdr3_candidates = list()
	translations = Seq(self.top_V.query_sequence)
	stops = bio.find_all_by_frames(translations, '*')
	for start in query_cdr3_starts:
		for end in query_cdr3_ends:
			log.verbose(f"{self.top_V.query_sequence[start:end+1]}, {start}, {end}: {start % 3 == (end-2) % 3}" , indent=2)
			candidate = CDR3(Seq(self.top_V.query_sequence[start:end+1]), start, end)
			candidate.get_translation()
			count_stops = len(stops[start % 3])
			candidate.calculate_score(count_stops)
			cdr3_candidates.append(candidate)
			log.verbose(f"{candidate.seq_aa}, {candidate.is_productive}, {candidate.score}", indent=3)
	return cdr3_candidates

def get_cdr3_2(self):
	start_query = self.top_V.reference_start - self.top_V.query_alignment_start
	abstract_V = self.top_V_reference.seq_nt[:start_query]
	abstract_J = self.top_J_reference.seq_nt[self.top_J.reference_end:]
	self.transcript = abstract_V + Seq(self.top_V.query_sequence) + abstract_J
	self.transcript_translations = bio.get_translation_frames(self.transcript)

	end_V = (start_query + self.top_V.query_alignment_end) // 3
	start_J = (start_query + self.top_J.query_alignment_start) // 3

	frame_stops = bio.find_all_by_frames(self.transcript_translations, '*')
	log.verbose(frame_stops)
	log.verbose(f"Checking for C, {end_V-12} (end_V-12) to {start_J}")
	cdr3_starts = bio.find_all(self.transcript_translations, 'C', end_V - 8, start_J)
	log.verbose(f"Checking for F, {end_V} to {start_J+12} (start_J+12)")
	cdr3_ends = bio.find_all(self.transcript_translations, 'F', end_V, start_J + 12)
	
	if not cdr3_starts:
		cdr3_starts = [end_V * 3 - 12]
	if not cdr3_ends:
		cdr3_ends = [start_J * 3 + 36]
	cdr3_candidates = list()
	for start in cdr3_starts:
		for end in cdr3_ends:
			candidate = CDR3(self.transcript[start:end+3], start, end+2)
			candidate.get_translation()
			if candidate.seq_aa is not None:
				count_stops = len(frame_stops[start % 3])
				candidate.calculate_score(count_stops)
			else:
				candidate.score = -999
			cdr3_candidates.append(candidate)
			
	self.cdr3 = max(cdr3_candidates, key=lambda c: c.score)
	log.verbose(f"WINNER: {self.cdr3.seq_aa}, {self.cdr3.is_productive}, {self.cdr3.score}", indent=2)
	print(self.VJ_alignment_diagram_full())

def get_cdr3_old1(self):
	log.info(self.query_name)
	log.info(self.VJ_alignment_stats())
	
	start_query = self.top_V.reference_start - self.top_V.query_alignment_start
	abstract_V = self.top_V_reference.seq_nt[:start_query]
	abstract_J = self.top_J_reference.seq_nt[self.top_J.reference_end:]
	self.transcript = abstract_V + Seq(self.top_V.query_sequence) + abstract_J
	self.transcript_translations = bio.get_translation_frames(self.transcript)

	cdr3_candidates = list()
	if self.top_V_reference.cdr3_positions is None \
		or bio.contains_indels(self.top_V) \
		or len(self.top_V_reference.seq_nt) - self.top_V.reference_end >= 12:
		query_cdr3_starts = bio.find_cdr3_start_via_aa(
			Seq(self.top_V.query_sequence[:(self.top_J.query_alignment_start-10)]), distance_from_end=40
		)
		if not query_cdr3_starts:
			query_cdr3_starts = [self.top_V.query_alignment_end - 20]
		log.verbose(f"Adjusted query_cdr3_starts! {query_cdr3_starts}")
	else:
		query_cdr3_starts = [
			self.top_V.query_alignment_start - self.top_V.reference_start + start \
				for start in self.top_V_reference.cdr3_positions
		]
	if self.top_J_reference.cdr3_positions is None \
		or bio.contains_indels(self.top_J):
		query_cdr3_ends = bio.find_cdr3_end_via_aa(
			Seq(self.top_V.query_sequence), index_start=self.top_V.query_alignment_end+10, length=40
		)
		if not query_cdr3_ends:
			query_cdr3_ends = [self.top_J.query_alignment_start + 40]
		log.verbose(f"Adjusted query_cdr3_ends! {query_cdr3_ends}")
	else:
		query_cdr3_ends = [
			self.top_J.query_alignment_start - self.top_J.reference_start + end \
				for end in self.top_J_reference.cdr3_positions
		]
	cdr3_candidates = self.score_cdr3_candidates(query_cdr3_starts, query_cdr3_ends)
	self.cdr3 = max(cdr3_candidates, key=lambda c: c.score)
	log.verbose(f"WINNER: {self.cdr3.seq_aa}, {self.cdr3.is_productive}, {self.cdr3.score}", indent=2)
	print(self.VJ_alignment_diagram_full())

def VJ_alignment_diagram_OLD(self) -> str:
	"""Includes just the part of the reference that aligns to the query"""
	return textwrap.dedent(
		f"""
		{' ' * self.top_V.query_alignment_start}{self.ref_seq_V[self.top_V.reference_start:self.top_V.reference_end]}
		{self.visualize_cigar(self.top_V.cigarstring)}
		{' ' * self.top_V.query_alignment_start}{'?'*(self.top_V.reference_end-self.top_V.reference_start)}
		{' ' * self.query_cdr3_start}*
		{self.top_V.query_sequence}
		{' ' * self.query_cdr3_end}*
		{' ' * self.top_J.query_alignment_start}{'?'*(self.top_J.reference_end-self.top_J.reference_start)}
		{self.visualize_cigar(self.top_J.cigarstring)}
		{' ' * self.top_J.query_alignment_start}{self.ref_seq_J[self.top_J.reference_start:self.top_J.reference_end]}
		"""
	)

def get_cdr3_positions_OLD(self, cdr3_positions: Dict[str, int]):
	self.ref_cdr3_start = cdr3_positions[self.top_V.reference_name] - 1 # Make 0 indexed
	self.query_cdr3_start = self.top_V.query_alignment_start - self.top_V.reference_start + self.ref_cdr3_start
	self.ref_cdr3_end = cdr3_positions[self.top_J.reference_name] - 1 # Make 0 indexed
	self.query_cdr3_end = self.top_J.query_alignment_start - self.top_J.reference_start + self.ref_cdr3_end 	
	# This code works well! It just doesn't really make much of a difference in the big picture.
	# I also don't know if it's worth potentially biasing against nDeletionsBefore%3!=0 
	# so I'm disabling for now.
	# Also just occurred to me that there may be deletions upstream of the query_sequence
	# sooo maybe reads with nDels % 3 != 0 can still be viable. Can't tell without ORF context.
	# I think it's worth talking to Sarah about what should be done here and for insertions.
	# Probably should just scrap top VJ and for the de Bruijn graph sequence assembly method anyways.
	"""
	if 'D' in self.top_V.cigarstring:
		deletion_positions = self.get_deletion_positions(self.top_V.cigarstring)
		deletions_before = 0
		for gap_start, gap_width in reversed(deletion_positions.items()):
			gap_end = gap_start + gap_width-1
			if gap_start <= self.query_cdr3_start <= gap_end:
				break
			elif gap_end < self.query_cdr3_start:
				deletions_before += gap_width
		if deletions_before % 3 == 0:
			self.query_cdr3_start -= deletions_before
	"""
	# I don't think adjusting the J end is as important. Might be hurtful for analysis.
	# This code needs rewriting if we want something like it.
	""" 
	if 'D' in self.top_J.cigarstring:
		deletion_positions = self.get_deletion_positions(self.top_J.cigarstring)
		for gap_start, gap_width in reversed(deletion_positions.items()):
			gap_end = gap_start + gap_width-1
			if gap_start <= self.query_cdr3_end <= gap_end:
				break
			if gap_end < self.query_cdr3_end - 2:
				self.query_cdr3_end -= gap_width
	"""
# DEPRECATED
def get_cdr3_sequence_OLD(self):
	if self.top_J.query_alignment_end <= self.top_V.query_alignment_end:
		#log.verbose(self.VJ_alignment_stats())
		#log.verbose(self.VJ_alignment_diagram(), indent=0)
		log.warn(
			f"End of J alignment ({self.top_J.query_alignment_end}) does not "
			f"come after the end of V alignment ({self.top_V.query_alignment_end})"
		)
	elif self.top_J.query_alignment_end >= self.query_cdr3_end: # If we have the entire J region
		self.is_complete_cdr3 = True
		self.cdr3 = self.top_V.query_sequence[self.query_cdr3_start:self.query_cdr3_end+1]
	else: # Incomplete CDR3
		self.cdr3 = self.top_V.query_sequence[self.query_cdr3_start:]
	"""#DEBUG
	for char in self.top_V.cigarstring:
		if 47 >= ord(char) or ord(char) >= 58:
			if char not in "SMDI":
				log.verbose(self.top_V.cigarstring)
	for char in self.top_J.cigarstring:
		if 47 >= ord(char) or ord(char) >= 58:
			if char not in "SMDI":
				log.verbose(self.top_J.cigarstring)
	"""

@log.fdoc
def VJ_alignment_stats_old(self) -> str:
	#TODO: Remove if new method succeeds
	return \
		f"""
		{'_'*80}
		| RNAME |{'{:^34}'.format(self.top_V.reference_name)}||{'{:^34}'.format(self.top_J.reference_name)}|
		| RBEG  |{'{:^34}'.format(self.top_V.reference_start)}||{'{:^34}'.format(self.top_J.reference_start)}|
		| REND  |{'{:^34}'.format(self.top_V.reference_end)}||{'{:^34}'.format(self.top_J.reference_end)}|
		| RCDR3 |{'{:^34}'.format(self.ref_cdr3_start)}||{'{:^34}'.format(self.ref_cdr3_end)}|
		|{' '*78}|
		| QBEG  |{'{:^34}'.format(self.top_V.query_alignment_start)}||{'{:^34}'.format(self.top_J.query_alignment_start)}|
		| QEND  |{'{:^34}'.format(self.top_V.query_alignment_end)}||{'{:^34}'.format(self.top_J.query_alignment_end)}|
		| QALEN |{'{:^34}'.format(self.top_V.query_alignment_length)}||{'{:^34}'.format(self.top_J.query_alignment_length)}|
		| QCDR3 |{'{:^34}'.format(self.query_cdr3_start)}||{'{:^34}'.format(self.query_cdr3_end)}|
		|{' '*78}|
		| CIGAR |{'{:^34}'.format(self.top_V.cigarstring)}||{'{:^34}'.format(self.top_J.cigarstring)}|
		| SCORE |{'{:^34}'.format(self.top_V.get_tag('AS'))}||{'{:^34}'.format(self.top_J.get_tag('AS'))}|
		| EDIST |{'{:^34}'.format(self.top_V.get_tag('NM'))}||{'{:^34}'.format(self.top_J.get_tag('NM'))}|
		| PHRED |{'{:^34}'.format(self.top_V.get_tag('UQ'))}||{'{:^34}'.format(self.top_J.get_tag('UQ'))}|
		{'_'*80}
		"""

	# NOTE: bowtie2 cigar doesn't accurately convey locations of mismatches.
	def visualize_cigar(self, cigar: str) -> str:
		output = ""
		num = ""
		for char in cigar:
			if  47 < ord(char) < 58: # If number 0-9
				num += char
			else:
				if char == 'S':
					sym = ' '
				elif char == 'M':
					sym = '|'
					# If equal then '|' else ' '
				elif char == 'D':
					sym = 'x'
				else: # I = 'I'
					sym = char
				output += f"{sym * int(num)}"
				num = ""
		return output