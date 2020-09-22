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

class Read(object):
	def __init__(self, query_name: str):
		self.query_name = query_name
		self.unique_subregions = set()
		self.alignments = list()
		self.top_J = None
		self.top_V = None
		#self.top_index_V = -1
		#self.top_index_J = -1
		self.is_TRA = None
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

	def __len__(self) -> int:
		return len(self.alignments)

	def __iter__(self) -> Iterator:
		return iter(self.alignments)

	def append(self, alignment: AlignedSegment):
		# If the alignment score is equal to the best alignment score for the Read,
		# bookmark that alignment as the Read's best alignment
		tags = dict(alignment.tags)
		if tags["AS"] == tags["XS"]:
			self.best_alignment = alignment
		self.alignments.append(alignment)

	@log.fdoc
	def print_alignment(self, alignment: AlignedSegment):
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
					score = alignment.get_tag("AS")
					if score >= top_score_V:
						if score == top_score_V:
							self.top_V = self.tie_break(self.top_V, alignment)
						else:
							self.top_V = alignment
							top_score_V = score
					break
			else:
				for j in ("TRAJ", "TRBJ"):
					if j in alignment.reference_name:
						self.has_region[j] = True
						score = alignment.get_tag("AS")
						if score >= top_score_J:
							if score == top_score_J:
								self.top_J = self.tie_break(self.top_J, alignment)
							else:
								self.top_J = alignment
								top_score_J = score
						break
				else:
					for c in ("TRAC", "TRBC"):
						if c in alignment.reference_name:
							self.has_region[c] = True
							break
					else:
						self.has_region["UNKN"] = True
		if "TRA" in self.top_V.reference_name and "TRA" in self.top_J.reference_name:
			self.is_TRA = True
		elif "TRB" in self.top_V.reference_name and "TRB" in self.top_J.reference_name:
			self.is_TRA = False
		else: 
			log.error("We have a top TRA/TRB combo!") #TODO: This is debug. If ever happens, what do?

	def get_cdr3_positions(self, cdr3_positions: Dict[str, int]):
		self.ref_cdr3_start = cdr3_positions[self.top_V.reference_name]
		self.query_cdr3_start = self.top_V.query_alignment_start + self.ref_cdr3_start - self.top_V.reference_start
		self.ref_cdr3_end = cdr3_positions[self.top_J.reference_name]
		self.query_cdr3_end = self.top_J.query_alignment_start + self.ref_cdr3_end - self.top_J.reference_start

	def get_cdr3_sequence(self):
		if self.top_J.query_alignment_end <= self.top_V.query_alignment_end:
			log.verbose(self.VJ_alignment_stats())
			log.verbose(self.VJ_alignment_diagram(), indent=0)
			log.warn(
				f"End of J alignment ({self.top_J.query_alignment_end}) does not "
				f"come after the end of V alignment ({self.top_V.query_alignment_end})"
			)
			return
		if self.top_J.query_alignment_end >= self.query_cdr3_end: # If we have the entire J region
			self.is_complete_cdr3 = True
			self.cdr3 = self.top_V.query_sequence[self.query_cdr3_start:self.query_cdr3_end]
		else: # Incomplete CDR3
			self.cdr3 = self.top_V.query_sequence[self.query_cdr3_start:]

	def VJ_alignment_diagram(self) -> str:
		return textwrap.dedent(
			f"""
			{' '*self.top_V.query_alignment_start}{self.ref_seq_V[self.top_V.reference_start:self.top_V.reference_end]}
			{' '*self.top_V.query_alignment_start}{'|'*(self.top_V.reference_end-self.top_V.reference_start)}
			{' '*self.query_cdr3_start}*
			{self.top_V.query_sequence}
			{' '*self.query_cdr3_end}*
			{' '*self.top_J.query_alignment_start}{'|'*(self.top_J.reference_end-self.top_J.reference_start)}
			{' '*self.top_J.query_alignment_start}{self.ref_seq_J[self.top_J.reference_start:self.top_J.reference_end]}
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
	
	#########################################################################
	#	Below is deprecated code that will be removed eventually
	#########################################################################

	def get_top_subregion_alignment(self, subregion: str) -> AlignedSegment:
		"""For use with _OLD_parse_bam"""
		top_score = 0
		top_alignment = None
		for alignment in self.alignments:
			if subregion in alignment.reference_name:
				self.unique_subregions.add(alignment.reference_name)
				tags = dict(alignment.tags)
				score = tags["AS"]
				best_score = tags["XS"]
				if score == best_score:
					log.verbose("Found best alignment for query")
					return alignment
				log.verbose(f"score for {alignment.reference_name} is {score}, will{' NOT' if score <= top_score else ''} beat {top_score}")	
				if score > top_score: # TODO: What happens if it ties?
					top_alignment = alignment
					top_score = score
					log.verbose(f"Beat score, now top region is {top_alignment.reference_name} with AS {top_score}")			
		return top_alignment