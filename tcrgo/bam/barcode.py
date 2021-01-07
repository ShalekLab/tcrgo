import pysam
import textwrap
from collections import Counter, defaultdict
from typing import List, Dict, Iterator, Set, Tuple
from .umi import UMI
from ..log import Log

log = Log("root")
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

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
