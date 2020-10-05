import pysam
import textwrap
from collections import Counter, defaultdict

# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set, Tuple, Counter
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile
from .read import Read

from ..io import Log
log = Log(name=__name__)
log.proceed()

class UMI(object):
	def __init__(self, sequence: str):
		self.sequence = sequence
		self.reads = dict()
		self.counts_region = Counter()
		#self.counts_top_V = Counter()
		#self.counts_top_J = Counter()
		#self.top_V = None
		#self.top_J = None
		self.top_VJ = None
		self.count_top_VJ = None
		self.frequency_top_VJ = None
		self.reads_top_VJ = None
		self.top_cdr3 = None
		self.is_complete_cdr3 = False
		self.count_top_cdr3 = None
		self.frequency_top_cdr3 = None

	def __len__(self) -> int:
		return len(self.reads)

	def __getitem__(self, key: str) -> Read:
		return self.reads[key]

	def __setitem__(self, key, value):
		self.reads[key] = value

	def __delitem__(self, key):
		del self.reads[key]

	def __iter__(self) -> Iterator:
		return iter(self.reads)

	def __contains__(self, item) -> bool:
		return item in self.reads

	def __missing__(self, key):
		raise KeyError(f"No Read {key} in UMI")

	def _keys(self) -> List[str]:
		return self.reads.keys()

	def _values(self) -> List[Read]:
		return self.reads.values()

	def items(self) -> List[Tuple[str, Read]]:
		return self.reads.items()

	def get_read_query_names(self) -> List[str]:
		return self._keys()

	def get_reads(self) -> List[Read]:
		return self._values() 

	def cdr3_statistics(self):
		count_complete = 0
		count_incomplete = 0
		count_total = 0
		log.info(f"UMI {self.sequence}, gathering CDR3 statistics")
		for read in self.get_reads():
			if read.cdr3_status == "COM":
				count_complete += 1
			elif read.cdr3_status == "INC":
				count_incomplete += 1
			count_total += 1
			#log.verbose(read.cdr3)
		log.verbose(
			f"For {self.sequence} got {count_complete} complete and "
			f"{count_incomplete} incomplete of {count_total} reads."
		)

	def find_top_VJ(self):
		counts_VJ = Counter()
		VJ_reads = defaultdict(list)
		for read in self.get_reads():
			VJ = f"{read.top_V.reference_name}|{read.top_J.reference_name}"
			VJ_reads[VJ].append(read)
			counts_VJ[VJ] += 1
		self.top_VJ = max(counts_VJ, key=counts_VJ.get)
		self.reads_top_VJ = VJ_reads[self.top_VJ]
		#self.top_V, self.top_J = self.top_VJ.split('|')

		self.count_top_VJ = counts_VJ[self.top_VJ]
		count_total = len(self)
		self.frequency_top_VJ = self.count_top_VJ / count_total
	
		#log.info(f"{counts_VJ=}")
		#log.verbose(f"{top_VJ}/{count_total} = {self.frequency_top_VJ}")
		#log.info(f"{count_total=}")

	def count_regions(self):
		for read in self.get_reads():
			for subregion in ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC", "UNKN"):
				if read.has_region[subregion]:
					self.counts_region[subregion] += 1

	def resolve_tcr_identity(self):
		"""
		Figures out the V, J, and CDR3 consensus sequences for the UMI
		TODO: Currently only considers NT, not AA sequence.
		TODO: Unfinished
		"""
		for read in self.reads_top_VJ:
			if read.cdr3 is not None and read.is_complete_cdr3:
				self.is_complete_cdr3 = True # UMI's CDR3 identity will only consider complete CDR3's
				break
		counts_cdr3 = Counter()
		for read in self.reads_top_VJ:
			if self.is_complete_cdr3 and not read.is_complete_cdr3:
				continue
			elif read.cdr3:
				counts_cdr3[read.cdr3] += 1
		if not counts_cdr3:
			log.warn(f"{self.sequence} had all CDR3s=None")
			self.top_cdr3 = "NA"
			self.count_top_cdr3 = "NA"
			self.frequency_top_cdr3 = "NA"
			return

		top_cdr3s = []
		max_count = 0
		for cdr3, count in counts_cdr3.items():
			if not top_cdr3s:
				top_cdr3s.append(cdr3)
				max_count = count
			elif count > max_count:
				top_cdr3s = [cdr3]
				max_count = count
			elif count == max_count:
				top_cdr3s.append(cdr3)
		if len(top_cdr3s) > 1:
			log.warn(f"We have a CDR3 tie of count {max_count} among {','.join(top_cdr3s)}")
			#TODO: Resolve ties via consensus sequence?
		top_cdr3 = top_cdr3s[0] # TODO: For now just do first appended. Resolve.
		
		''' TODO: Hamming variants
		cdr3_candidates = []
		for cdr3 in counts_cdr3.keys():
			if len(cdr3) == len(top_cdr3):
		'''
		self.count_top_cdr3 = counts_cdr3[top_cdr3]
		self.frequency_top_cdr3 = self.count_top_cdr3 / len(self)
		if not self.is_complete_cdr3:
			top_cdr3 += '_'
		self.top_cdr3 = top_cdr3

	#########################################################################
	#	Below is deprecated code that will be removed eventually
	#########################################################################

	""" # This was a more general approach
	def count_subregions(self):
		for read in self.get_reads():
			for subregion in read.unique_subregions:
				if "TRAV" in subregion:
					self.subregion_counts["TRAV"] += 1
				elif "TRAJ" in subregion:
					self.subregion_counts["TRAV"] += 1
				elif "TRBV" in subregion:
					self.subregion_counts["TRAV"] += 1
				elif "TRBJ" in subregion:
					self.subregion_counts["TRAV"] += 1
				elif "TRAC" in subregion:
					self.subregion_counts["TRAV"] += 1
				elif "TRBC" in subregion:
					self.subregion_counts["TRAV"] += 1
				else:
					self.subregion_counts["TRAV"] += 1

	# This was a more specific approach
	def count_subregions(self):
		for read in self.get_reads():
			for subregion in read.unique_subregions:
				if 'AV' in subregion:
					self.count_AV[subregion] += 1
				elif 'AJ' in subregion:
					self.count_AJ[subregion] += 1
				elif 'BV' in subregion:
					self.count_BV[subregion] += 1
				elif 'BJ' in subregion:
					self.count_BJ[subregion] += 1
	"""

	def calculate_subregion_statistics(self) -> str:
		sum_AV = sum(self.count_AV.values())
		sum_AJ = sum(self.count_AJ.values())
		sum_BV = sum(self.count_BV.values())
		sum_BJ = sum(self.count_BJ.values())
		total = sum_AV + sum_AJ + sum_BV + sum_BJ
		return f"{sum_AV=}, {sum_AV/total}; {sum_AJ=}, {sum_AJ/total};" \
				f"{sum_BV=}, {sum_BV/total}; {sum_BJ=}, {sum_BJ/total}"

	def discard_single_alignments(self):
		"""DEPRECATED."""
		to_delete = []
		for read in self:
			if len(read) == 1:
				log.verbose(f"Deleting {read.query_name}")
				to_delete.append(read)
		for read in to_delete:
			del read

	@log.time
	def find_top_V_J_separately(self) -> Tuple[Counter[str], Counter[str]]:
		counts_top_V = Counter()
		counts_top_J = Counter()
		for read in self.get_reads():
			counts_top_V[read.top_V.reference_name] += 1
			counts_top_J[read.top_J.reference_name] += 1
		top_V = max(counts_top_V, key=counts_top_V.get)
		top_J = max(counts_top_J, key=counts_top_J.get)
		count_total = sum(counts_top_V.values())
		frequency_top_V = counts_top_V[top_V] / count_total
		frequency_top_J = counts_top_J[top_J] / count_total
		log.info(f"{counts_top_V=}")
		log.verbose(f"{top_V=} {frequency_top_V=}")
		log.info(f"{counts_top_J=}")
		log.verbose(f"{top_J=} {frequency_top_J=}")
		log.info(f"{count_total=}")
		return 

	@log.time
	def find_top_V_J_by_alignments(self):
		self.counts_top_V = Counter()
		self.counts_top_J = Counter()
		count_total_J = 0
		count_total_V = 0
		for read in self.get_reads():
			for subregion in read.unique_subregions:
				if "TRAV" in subregion or "TRBV" in subregion:
					self.counts_top_V[subregion] += 1
					count_total_V += 1
				elif "TRAJ" in subregion or "TRBJ" in subregion:
					self.counts_top_J[subregion] += 1
					count_total_J += 1
				
		top_V = max(self.counts_top_V, key=self.counts_top_V.get)
		top_J = max(self.counts_top_J, key=self.counts_top_J.get)
		
		frequency_top_V = self.counts_top_V[top_V] / count_total_V
		frequency_top_J = self.counts_top_J[top_J] / count_total_J
		
		log.info(f"{self.counts_top_V=}")
		log.verbose(f"{top_V=} {frequency_top_V=}")
		log.info(f"{self.counts_top_J=}")
		log.verbose(f"{top_J=} {frequency_top_J=}")
		log.info(f"{count_total_V=} {count_total_J=}")