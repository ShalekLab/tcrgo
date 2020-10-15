import pysam
import textwrap

# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from collections import Counter, defaultdict
from typing import List, Dict, Iterator, Set, Tuple
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile
from .umi import UMI
from ..bio import levenshtein_distance

from ..log import Log
log = Log(name=__name__)
log.proceed()

class Barcode(object):
	def __init__(self, sequence: str):
		self.sequence = sequence
		self.umis = {}

	def __len__(self) -> int:
		return len(self.umis)

	def __getitem__(self, key: str) -> UMI:
		return self.umis[key]

	def __setitem__(self, key, value):
		self.umis[key] = value

	def __delitem__(self, key):
		del self.umis[key]

	def __iter__(self) -> Iterator:
		return iter(self.umis)

	def __contains__(self, item) -> bool:
		return item in self.umis

	def __missing__(self, key):
		raise KeyError(f"No UMI {key} in Barcode")

	def _keys(self) -> List[str]:
		return self.umis.keys()

	def _values(self) -> List[UMI]:
		return self.umis.values()

	def items(self) -> List[Tuple[str, UMI]]:
		return self.umis.items()

	def get_umi_sequences(self) -> List[str]:
		return self._keys()

	def get_umis(self) -> List[UMI]:
		return self._values()

	def resolve_cdr3s(self):
		#for other_chain in ("TRA", "TRB"):
		leaders = Counter()
		candidates = Counter() # In case no leaders
		umis_contested = list()
		cdr3s_contested = list()
		for umi in self.get_umis():
			if umi.top_cdr3 is None: #or other_chain in umi.top_VJ:
				continue
			if len(umi.top_cdr3) == 0:
				umi.top_cdr3 = None
			elif len(umi.top_cdr3) == 1:
				umi.top_cdr3 = umi.top_cdr3.pop()
				leaders[umi.top_cdr3] += 1
				candidates[umi.top_cdr3] += 1
			else:
				umis_contested.append(umi)
				cdr3s_contested.append(umi.top_cdr3)
				for cdr3 in umi.top_cdr3:
					candidates[cdr3] += 1
		return # DEBUG.
		"""
		Giving up on this for now
		"""
		if not leaders:
			for umi in self.get_umis():
				log.verbose(umi.top_cdr3, indent=2)
			log.verbose(leaders)
			log.verbose(candidates)
			log.verbose(umis_contested)
			log.verbose(cdr3s_contested)
			log.error("Well it happened")
			#leaders = candidates

		#cdr3_first, cdr3_second = leaders.most_common(2)
		leaders_sorted = [k for k,v in sorted(leaders, key=lambda item: item[1], reverse=True)]
		for i in range(len(umis_contested)):
			count_best = 0
			cdr3_best = None
			for cdr3 in cdr3s_contested[i]:
				for cdr3_leader in leaders_sorted:
					if levenshtein_distance(cdr3_leader, cdr3) <= 1:
						leader_count = leaders[cdr3_leader]
						if leader_count > count_best:
							count_best = leader_count
							cdr3_best = cdr3
						break
			if cdr3_best != None:
				umis_contested[i].cdr3 = cdr3_best
			else:
				log.verbose(leaders)
				log.verbose(candidates)
				log.verbose(umis_contested)
				log.verbose(cdr3s_contested)
				log.error("Still couldn't resolve")
		

		

