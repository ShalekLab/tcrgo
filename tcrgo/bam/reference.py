from typing import Tuple, List, Dict, Iterator, Union, Any, Set
from Bio.Seq import Seq
from itertools import chain, combinations
import tcrgo.bio as bio
from ..log import Log
from ..io import parse_fasta

log = Log("root")

class Reference(object):
	def __init__(self, name: str, seq_nt: Union[str, Seq], 
		cdr3_positions: List[int]=None, orf: int=None):
		self.name = name
		self.root = self.name.split('-')[0]
		self.orf = orf # None for C + J
		if isinstance(seq_nt, str):
			self.seq_nt = Seq(seq_nt)
		else:
			self.seq_nt = seq_nt
		self.seq_aa = None
		# I should have made subclasses for V, J, and C, or made all List[int]
		# would get rid of the problem below and tidy up a lot of this module.
		self.cdr3_positions = cdr3_positions # int for 'V', ((int?),(int?),(int?)) for 'J'
		self.is_TRA, self.segment = self.get_identity(name)

	def get_identity(self, name: str) -> Tuple[bool, str]:
		is_TRA = False
		if name.startswith("TRA"):
			is_TRA = True
		elif not name.startswith("TRB"):
			log.error(f"{name} does not begin with 'TRA' or 'TRB'!")
		for segment in ('V', 'J', 'C'):
			if name[3] == segment:
				return is_TRA, name[3]
		else:
			log.error(f"{name} could not be identified as 'V', 'J', or 'C' segment "
				f"(got {name[3]} as 4th character of reference name)")

	def has_cdr1(self, seqs_aa: Tuple[Seq, ...]) -> List[bool]:
		cdr1_starts = bio.find_all_in_frames(seqs_aa, 'C', 19, 24)
		cdr1_ends = bio.find_all_in_frames(seqs_aa, 'W', 28, 44)
		has_cdr1 = [False, False, False]
		for f in range(3):
			if cdr1_starts[f] and cdr1_ends[f]:
				has_cdr1[f] = True
		return has_cdr1

	def identify_best_frame(self):
		"""If no cdr3 positions file, automatically determine"""
		if self.segment == 'C':
			return
		seqs_aa = bio.get_frame_translations(self.seq_nt)
		if self.segment == 'J':
			self.cdr3_positions = bio.find_all_in_frames(seqs_aa, 'F', 5, 12)
			return
		frames_has_cdr1 = self.has_cdr1(seqs_aa)
		count_has_cdr1 = frames_has_cdr1.count(True)
		starts = bio.find_all_in_frames(seqs_aa, 'C', 6, are_from_end=True)
		if count_has_cdr1 == 1:
			self.orf = frames_has_cdr1.index(True)
			if not starts[self.orf]:
				self.cdr3_positions = \
					len(self.seq_nt) - (len(self.seq_nt) - self.orf) % 3 - 15
			else:
				self.cdr3_positions = min(starts[self.orf])
			return
		has_start = [len(frame) > 0 for frame in starts]
		count_has_start = has_start.count(True)
		if count_has_start == 1:
			self.orf = has_start.index(True)
			self.cdr3_positions = min(starts[self.orf])
			return	
		stops = bio.find_all_in_frames(seqs_aa, '*')
		frames_best_stop = bio.get_frame_min(stops)
		if count_has_start == 0 :
			self.orf = frames_best_stop.pop()
			self.cdr3_positions = \
				len(self.seq_nt) - (len(self.seq_nt) - self.orf) % 3 - 15
			return	
		frames_best_start = bio.get_frame_max(starts)
		frames_best = frames_best_stop & frames_best_start
		if not frames_best:
			frames_best = frames_best_stop | frames_best_start
		self.orf = frames_best.pop()
		self.cdr3_positions = min(starts[self.orf])

	def _verify_cdr3_position(self, minimum: int, maximum: int, 
		boundary: str, is_zero_indexed: bool=True):
		is_inaccessible = False
		if self.cdr3_positions < minimum or maximum < self.cdr3_positions:
			index = "Zero-indexed" if is_zero_indexed else "One-indexed"
			log.warn(
				f"{index} CDR3 {boundary} position {self.cdr3_positions} is out of "
				f"range [{minimum},{maximum}] for {self.name} ({len(self.seq_nt)}nt)",
				indent=1
			)
			is_inaccessible = True
		if not is_zero_indexed:
			self.cdr3_positions -= 1
		self.orf = (self.cdr3_positions + 3) % 3
		if self.segment == 'J':
			frames_positions = [[], [], []]
			frames_positions[self.orf].append(self.cdr3_positions)
			self.cdr3_positions = tuple([tuple(frame) for frame in frames_positions])
		return is_inaccessible

	def _verify_cdr3_positions(self, minimum: int, maximum: int, 
		boundary: str, is_zero_indexed: bool=True):
		is_inaccessible = False
		for cdr3_position in self.cdr3_positions:
			if cdr3_position < minimum or maximum < cdr3_position:
				index = "Zero-indexed" if is_zero_indexed else "One-indexed"
				log.warn(
					f"{index} CDR3 {boundary} position {cdr3_position} is out of "
					f"range [{minimum},{maximum}] for {self.name} ({len(self.seq_nt)}nt)",
					indent=1
				)
				is_inaccessible = True
		if not is_zero_indexed:
			self.cdr3_positions = [cdr3_position - 1 for cdr3_position in self.cdr3_positions]
		return is_inaccessible

	def verify_cdr3_positions(self, is_zero_indexed: bool=True):
		if self.segment == 'V':
			minimum = 1
			maximum = len(self.seq_nt)-2
			boundary = "starting"
		elif self.segment == 'J':
			minimum = 3
			maximum = len(self.seq_nt)
			boundary = "ending"
		else:
			return
		if is_zero_indexed:
			minimum -= 1
			maximum -= 1
		if isinstance(self.cdr3_positions, int):
			is_inaccessible = self._verify_cdr3_position(minimum, maximum, boundary, is_zero_indexed)
		elif isinstance(self.cdr3_positions, list):
			is_inaccessible = self._verify_cdr3_positions(minimum, maximum, boundary, is_zero_indexed)
		else:
			log.error("self.cdr3_positions must be int or Tuple[int], "
				f"got type {type(self.cdr3_positions)}.")
		return is_inaccessible

class ReferenceDict(object):
	def __init__(self, fasta: str=None):
		self.references = dict()
		if fasta is not None:
			self.build(parse_fasta(fasta))

	def __len__(self) -> int:
		return len(self.references)

	def __getitem__(self, key: str) -> Reference:
		return self.references[key]

	def __setitem__(self, key, value):
		self.references[key] = value

	def __delitem__(self, key):
		del self.references[key]

	def __iter__(self) -> Iterator:
		return iter(self.references)

	def __contains__(self, item) -> bool:
		return item in self.references

	def __missing__(self, key):
		raise KeyError("No Reference in ReferenceDict")

	def keys(self) -> List[str]:
		return self.references.keys()

	def values(self) -> List[Reference]:
		return self.references.values()

	def items(self) -> List[Tuple[str, Reference]]:
		return self.references.items()

	def __str__(self) -> str:
		if not self.values():
			return "Empty"
		return '\n'.join([ \
			f"{ref.name} - len {len(ref.seq_nt)} nt, frame {ref.orf}, cdr3pos {ref.cdr3_positions}" \
				for ref in self.values()
		])

	def build(self, entries: List[Tuple[str, str]], \
		cdr3_positions: Dict[str, int]=None, \
		is_zero_indexed: bool=False) -> Dict[str, Reference]:
		"""Build a dictionary of References from entries and sequences"""
		references = dict()
		count_inaccessible = 0
		for entry, seq_nt in entries:
			position = None
			ref = Reference(entry, seq_nt, position)
			if ref.segment != 'C':
				if cdr3_positions is not None:
					ref.cdr3_positions = cdr3_positions[entry]
					if ref.verify_cdr3_positions(is_zero_indexed):
						count_inaccessible += 1
				else:
					ref.identify_best_frame()
			references[entry] = ref
		self.references = references
		if count_inaccessible > 0:
			log.error(f"{count_inaccessible} out of range positions, see above "
				"warnings! Use --zero-indexed argument if not one-indexed. "
				"Or run without --cdr3-positions-file to auto-determine.")
		return references
		
	def verify(self):
		"""Ensure no duplicate entries, warn of V segment ambiguity"""
		for ref1, ref2 in combinations(list(self.references.values()), 2):
			if ref1.name == ref2.name:
				log.error(f"{ref1.name} and {ref2.name} have same entry name!")
			if ref1.seq_nt == ref2.seq_nt:
				log.error(f"{ref1.name} and {ref2.name} have same sequence!")
			if ref1.segment == 'V' and ref2.segment == 'V' and ref1.is_TRA == ref2.is_TRA:
				try:
					distance = bio.hamming_distance(ref1.seq_nt[-90:], ref2.seq_nt[-90:])
					if distance <= 3:
						log.warn(f"Last 90 bases of {ref1.name} and {ref2.name}: " 
						f"Hamming Distance={distance}", indent=1)
				except ValueError: # Sequences with len < 90 will be skipped.
					pass