import pysam
import textwrap
from collections import Counter

# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set, Tuple
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
		self.count_AV = Counter()
		self.count_AJ = Counter()
		self.count_BV = Counter()
		self.count_BJ = Counter()
		#self.unique_cdr3 = set()

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

	def _items(self) -> List[Tuple[str, Read]]:
		return self.reads.items()

	def get_read_query_names(self) -> List[str]:
		return self._keys()

	def get_reads(self) -> List[Read]:
		return self._values()

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

	def calculate_subregion_statistics(self) -> str:
		sum_AV = sum(self.count_AV.values())
		sum_AJ = sum(self.count_AJ.values())
		sum_BV = sum(self.count_BV.values())
		sum_BJ = sum(self.count_BJ.values())
		total = sum_AV + sum_AJ + sum_BV + sum_BJ
		return f"{sum_AV=}, {sum_AV/total}; {sum_AJ=}, {sum_AJ/total};" \
				f"{sum_BV=}, {sum_BV/total}; {sum_BJ=}, {sum_BJ/total}"

	def resolve_tcr_identity(self):
		"""Figures out the V, J, and CDR3 consensus sequences"""
		pass

	def discard_single_alignments(self):
		"""DEPRECATED."""
		to_delete = []
		for read in self:
			if len(read) is 1:
				log.verbose(f"Deleting {read.query_name}")
				to_delete.append(read)
		for read in to_delete:
			del read

	