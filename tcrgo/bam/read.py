import pysam
import textwrap

# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set
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
		self.best_alignment = None
		self.top_J = None
		self.top_V = None
		self.is_TRA = None
		self.cdr3 = None

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

	# TODO: Report ties for top score and break ties based on
	#		mapping quality, if same then alphabetical name of region
	def get_alignments(self, bam_indexed: IndexedReads):
		top_score_J = 0
		top_score_V = 0
		count = 0
		for alignment in bam_indexed.find(self.query_name):
			self.alignments.append(alignment)
			#self.unique_subregions.add(alignment.reference_name)
			#log.verbose(f"{alignment.reference_name}")
			score = dict(alignment.tags)["AS"]
			if 'V' in alignment.reference_name:
				if score > top_score_V:
					self.top_V = alignment
					top_score_V = score
					#log.verbose(f"Top V is now {read.top_V.reference_name} with score={top_score_V}", level=2)
			elif 'J' in alignment.reference_name:
				if score > top_score_J:
					self.top_J = alignment
					top_score_J = score
					#log.verbose(f"Top J is now {read.top_J.reference_name} with score={top_score_J}", level=2)
			if 'A' in alignment.reference_name:
				self.is_TRA = True
			else: #TRB
				self.is_TRA = False

	def get_top_subregion_alignment(self, subregion: str) -> AlignedSegment:
		"""For use with parse_bam1"""
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

	def get_cdr3_position(self, alignment: AlignedSegment, cdr3_positions: Dict[str, int]):
		log.info(f"Getting CDR3 base pos for {self.query_name} {alignment.reference_name}.")
		return cdr3_positions[alignment.reference_name.replace("macfa", '')] # TODO find better way to account for fasta and CDR3 name differences

	"""
	def __str__(self) -> str:
		return "	READ" + textwrap.indent(
			textwrap.dedent(
				#f
				#QNAME:	{self.query_name}
				#NALIGN:	{len(self.alignments)}
				#SUBREG:	{self.unique_subregions}
				#TOPV:	{self.top_V.reference_name if self.top_V is not None else self.top_V}
				#TOPJ:	{self.top_J.reference_name if self.top_J is not None else self.top_J}
				#
			),
			prefix='	'*3
		)
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
	
