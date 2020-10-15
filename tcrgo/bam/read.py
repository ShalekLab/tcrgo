import pysam
import textwrap
import re
from collections import defaultdict, Counter

# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set, Tuple
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

from ..log import Log
log = Log(name=__name__)
log.proceed()

import os.path as osp # DEV

class Read(object):
	def __init__(self, query_name: str):
		self.query_name = query_name
		self.unique_subregions = set()
		self.alignments = list()
		self.top_J = None
		self.top_V = None
		#self.top_index_V = -1 # for explore.py
		#self.top_index_J = -1 # for explore.py
		self.cdr3 = None
		self.is_complete_cdr3 = False
		self.has_region = defaultdict(bool)
		self.ties = Counter()
		# TODO: Consider getting rid of these and finding
		# way to reconfigure VJ_alignment_stats+_diagram
		self.ref_cdr3_start = None
		self.query_cdr3_start = None
		self.ref_cdr3_end = None
		self.query_cdr3_end = None
		self.ref_seq_J = ""
		self.ref_seq_V = ""

		self.scores = dict()
		self.cigars = dict()
		self.edit_distances = dict()
		self.mapqs = dict()
		self.phreds = dict()

	def __len__(self) -> int:
		return len(self.alignments)

	def __iter__(self) -> Iterator:
		return iter(self.alignments)

	def append(self, alignment: AlignedSegment):
		self.alignments.append(alignment)

	@log.fdoc
	def alignment_info(self, alignment: AlignedSegment) -> str:
		if alignment == None:
			return "NONE!"
		return \
			f"""
			QNAME:	{self.query_name}
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

	# TODO: Refine tie-breaking conditions
	# TODO: Output the tiebreak info as a file
	def tie_break(self, top_alignment: AlignedSegment, alignment: AlignedSegment) -> AlignedSegment:
		if alignment.get_tag("NM") < top_alignment.get_tag("NM"): # Edit Distance
			winner = alignment
			loser = top_alignment
			method = "EditDistance"
		elif alignment.get_tag("UQ") > top_alignment.get_tag("UQ"): # Phred
			winner = alignment
			loser = top_alignment
			method = "Phred"
		else: # Alphabetical order of reference name
			winner, loser = sorted((top_alignment, alignment), key=lambda a: a.reference_name)
			method = "Alphabetical"
		self.ties[f"{winner}>{loser}:{method}"] += 1
		return winner

	def compare_scores(self, top_alignment: AlignedSegment, top_score: int, alignment: AlignedSegment) -> Tuple[AlignedSegment, int]:
		"""Compare alignment scores and return the best alignment and the new best score."""
		score = alignment.get_tag("AS")
		if score >= top_score:
			if score == top_score:
				top_alignment = self.tie_break(top_alignment, alignment)
			else:
				top_alignment = alignment
				top_score = score
		return top_alignment, top_score

	# TODO: Normalize AS by length of ??? ?
	def parse_alignments(self, bam_indexed: IndexedReads):
		top_score_J = 0
		top_score_V = 0
		for alignment in bam_indexed.find(self.query_name):
			self.alignments.append(alignment)
			self.unique_subregions.add(alignment.reference_name)
			for v in ("TRAV", "TRBV"):
				if v in alignment.reference_name:
					self.has_region[v] = True
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
		"""
		if "TRA" in self.top_V.reference_name and "TRA" in self.top_J.reference_name:
			
		elif "TRB" in self.top_V.reference_name and "TRB" in self.top_J.reference_name:
			
		else: #TODO: This is debug. If ever happens, what do?
		"""
			""" # DEV: For writing TRA/TRB combos to file
			if not osp.isfile("TRABcombos.tsv"):
				with open("TRABcombos.tsv", 'w') as combos:
					combos.write('\t'.join(["QNAME", "REF", "SCORE", "CIGAR", "PHRED", "EDIST", "FLAG", "MAPQ"])+'\n')
			with open("TRABcombos.tsv", 'a') as combos:
				combos.write(self.alignment_info_line(self.top_V))
				combos.write(self.alignment_info_line(self.top_J))
			"""

	def get_cdr3_positions(self, cdr3_positions: Dict[str, int]):
		""" TODO: I don't think this handles deletions! Check
		if 'D' in self.top_V.cigarstring:
			deletion_positions = [distinct deletion positions in query]
			deletion_gaps = [number of bases deleted]
			for i in range(len(deletion_positions)):
				if deletion_positions[i] + deletion_gaps[i]-1 <= cdr3_query_start:
					cdr3_query_start -= 1
					cdr3_query_end -= 1
				elif deletion_positions[i] + deletion_gaps[i]-1 <= cdr3_query_end:
					cdr3_query_end -= 1
		"""
		self.ref_cdr3_start = cdr3_positions[self.top_V.reference_name] - 1 # Make 0 indexed
		self.query_cdr3_start = self.top_V.query_alignment_start - self.top_V.reference_start + self.ref_cdr3_start
		self.ref_cdr3_end = cdr3_positions[self.top_J.reference_name] - 1 # Make 0 indexed
		self.query_cdr3_end = self.top_J.query_alignment_start - self.top_J.reference_start + self.ref_cdr3_end 

	def get_cdr3_sequence(self):
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
		#log.verbose(self.VJ_alignment_stats())
		#log.verbose(self.VJ_alignment_diagram_full(), indent=0)
		#log.verbose(self.cdr3)

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
				elif char == 'D':
					sym = 'x'
				else:
					sym = char
				output += f"{sym * int(num)}"
				num = ""
		return output

	def VJ_alignment_diagram(self) -> str:
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

	def VJ_alignment_diagram_full(self) -> str:
		query_start = self.top_V.reference_start - self.top_V.query_alignment_start
		j_start = query_start + self.top_J.query_alignment_start - self.top_J.reference_start
		return textwrap.dedent(
			f"""
			{self.ref_seq_V}
			{' ' * self.top_V.reference_start}{self.ref_seq_V[self.top_V.reference_start:self.top_V.reference_end]}
			{' ' * query_start}{self.visualize_cigar(self.top_V.cigarstring)}
			{' ' * self.ref_cdr3_start}r
			{' ' * query_start}{' ' * self.query_cdr3_start}q
			{' ' * query_start}{self.top_V.query_sequence}
			{' ' * query_start}{' ' * self.query_cdr3_end}q
			{' ' * j_start}{' '*self.ref_cdr3_end}r
			{' ' * query_start}{self.visualize_cigar(self.top_J.cigarstring)}
			{' ' * query_start}{' '*self.top_J.query_alignment_start}{self.ref_seq_J[self.top_J.reference_start:self.top_J.reference_end]}
			{' ' * j_start}{self.ref_seq_J}
			"""
		)

	@log.fdoc
	def VJ_alignment_stats(self) -> str:
		return \
			f"""
			{'_'*36}
			{'{:^36}'.format(self.top_V.reference_name)}
			{self.top_V.reference_start=}
			{self.top_V.reference_end=}
			{self.ref_cdr3_start=}
			-
			{self.top_V.query_alignment_start=}
			{self.top_V.query_alignment_end=}
			{self.query_cdr3_start=}
			-
			{self.top_V.cigarstring}
			{'_'*36}
			{'{:^36}'.format(self.top_J.reference_name)}
			{self.top_J.reference_start=}
			{self.top_J.reference_end=}
			{self.ref_cdr3_end=}
			-
			{self.top_J.query_alignment_start=}
			{self.top_J.query_alignment_end=}
			{self.query_cdr3_end=}
			-
			{self.top_J.cigarstring}
			{'_'*36}
			"""

	@log.fdoc
	def __str__(self) -> str:
		return \
			f"""
			QNAME:	{self.query_name}
			NALIGN:	{len(self.alignments)}
			SUBREG:	{self.unique_subregions}
			TOPV:	{self.top_V.reference_name if self.top_V is not None else self.top_V}
			TOPJ:	{self.top_J.reference_name if self.top_J is not None else self.top_J}
			"""
