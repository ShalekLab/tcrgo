from typing import Tuple, List, Dict, Iterator, Union, Any, Set
from Bio.Seq import Seq
from itertools import chain, combinations
import tcrgo.bio as bio
from ..log import Log

log = Log(name=__name__)
log.proceed()

class Reference(object):
	def __init__(self, name: str, seq_nt: Union[str, Seq], 
		cdr3_positions: List[int]=None, orf: int=None):
		self.name = name
		self.root = self.name.split('-')[0]
		self.orf = orf
		if isinstance(seq_nt, str):
			self.seq_nt = Seq(seq_nt)
		else:
			self.seq_nt = seq_nt
		self.seq_aa = None
		self.cdr3_positions = cdr3_positions
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


	def identify_best_frame(self):
		if self.segment == 'C':
			return
		seqs_aa = bio.get_frame_translations(self.seq_nt)
		if self.segment == 'J':
			self.cdr3_positions = bio.find_all_in_frames(seqs_aa, 'F', 5, 12)
			return
		has_cdr1 = bio.has_cdr1(seqs_aa)
		count_has_cdr1 = has_cdr1.count(True)
		starts = bio.find_all_in_frames(seqs_aa, 'C', 6, are_from_end=True)
		if count_has_cdr1 == 1:
			self.orf = has_cdr1.index(True)
			if not starts[self.orf]:
				self.cdr3_positions = \
					len(self.seq_nt) - (len(self.seq_nt) - self.orf) % 3 - 15
			else:
				self.cdr3_positions = min(starts[self.orf])
			return	
		has_start = [len(frame) > 0 for frame in starts]
		count_has_start = has_start.count(True)
		print(count_has_start)
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
		print(self.name, self.orf, self.cdr3_positions)

	def find_cdr3_start_via_aa(self, distance_from_end: int=6):
		aa_sequences = dict()
		for frame in (0,1,2):
			end = len(self.seq_nt) - ((len(self.seq_nt) - frame) % 3)
			aa_sequences[frame] = self.seq_nt[frame:end].translate()
		positions = list()
		for frame, seq_nt in aa_sequences.items():
			if seq_nt.find('*') == -1:
				start = len(seq_nt) - distance_from_end
				end = len(seq_nt) # endslice
				while start < end:
					position = seq_nt.find('C', start, end)
					if position == -1:
						break
					positions.append(position * 3 + frame)	
					start = position + 1
		if not positions:
			log.warn(f"No starting CDR3 cysteine codon detected for {self.name}", indent=1)
			for frame, seq_nt in aa_sequences.items():
				print(seq_nt, frame)
		else:
			self.cdr3_positions = sorted(positions, reverse=False)

	def find_cdr3_start(self, distance_end: int=0, distance_start: int=18):
		"""Greedily searches for TGT and TGC. The window to search in is emperical."""
		positions = {"TGT": list(), "TGC": list()}
		for codon in positions.keys():
			end = len(self.seq_nt) - distance_end # No -1 since end slice
			start = len(self.seq_nt) - distance_start - 1
			while end >= start + 3:
				position = self.seq_nt.rfind(codon, start, end)
				if position == -1:
					break
				positions[codon].append(position)	
				end = position + 2 # Because endslice
		if not positions["TGT"] and not positions["TGC"]:
			log.warn(f"No starting CDR3 codon detected for {self.name}", indent=1)
		else:
			self.cdr3_positions = sorted(positions["TGT"] + positions["TGC"], reverse=False)

	def find_cdr3_end_via_aa(self, index_start=4, distance_from_start: int=15):
		aa_sequences = dict()
		for frame in (0,1,2):
			end = len(self.seq_nt) - ((len(self.seq_nt) - frame) % 3)
			aa_sequences[frame] = self.seq_nt[frame:end].translate()
		positions = list()
		for frame, seq_nt in aa_sequences.items():
			if seq_nt.find('*') == -1:
				start = index_start
				end = index_start + distance_from_start # endslice
				while start < end:
					position = seq_nt.find('F', start, end)
					if position == -1:
						break
					positions.append(position * 3 + frame + 2)	
					start = position + 1
		if not positions:
			log.warn(f"No ending CDR3 phenylalanine codon detected for {self.name}", indent=1)
			#for frame, seq_nt in aa_sequences.items():
			#	print(seq_nt, frame)
		else:
			self.cdr3_positions = sorted(positions, reverse=True)

	def find_cdr3_end(self, index_start: int=14, length: int=23):
		"""Greedily searches for TTT and TTC. The window to search in is emperical."""
		positions = {"TTT": list(), "TTC": list()}
		for codon in positions.keys():
			start = index_start
			end = index_start + length # -1 + 1 since endslice
			while start <= end - 3:
				position = self.seq_nt.find(codon, start, end)
				if position == -1:
					break
				positions[codon].append(position + 2)
				start = position + 1
		if not positions["TTT"] and not positions["TTC"]:
			log.warn(f"No ending CDR3 codon detected for {self.name}", indent=1)
		else:
			self.cdr3_positions = sorted(positions["TTT"] + positions["TTC"], reverse=True)
			
	def find_cdr3_position(self):
		if self.segment == 'V':
			self.find_cdr3_start_via_aa()
		elif self.segment == 'J':
			self.find_cdr3_end_via_aa()
		else:
			log.warn(f"{self.name} is not a V or J segment.")

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
		for cdr3_position in self.cdr3_positions:
			if cdr3_position < minimum or maximum < cdr3_position:
				index = "Zero-indexed" if is_zero_indexed else "One-indexed"
				log.warn(
					f"{index} CDR3 {boundary} position {cdr3_position} is out of "
					f"range [{minimum},{maximum}] for {self.name} ({len(self.seq_nt)}nt)",
					indent=1
				)
		if not is_zero_indexed:
			self.cdr3_positions = [cdr3_position-1 for cdr3_position in self.cdr3_positions]

class ReferenceDict(object):
	def __init__(self):
		self.references = dict()
	
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

	def build(self, entries: List[Tuple[str, str]], \
		cdr3_positions: Dict[str, int]=None, \
		is_zero_indexed: bool=False) -> Dict[str, Reference]:
		"""Build a dictionary of References from entries and sequences"""
		for entry, seq_nt in entries:
			position = None
			ref = Reference(entry, seq_nt, position)
			if ref.segment != 'C':
				if cdr3_positions is not None:
					ref.cdr3_position = cdr3_positions[entry]
					ref.verify_cdr3_positions(is_zero_indexed)
				else:
					#ref.find_cdr3_position()
					ref.identify_best_frame()
			self.references[entry] = ref
	
	def verify(self):
		"""Ensure no duplicate entries, warn of V segment ambiguity"""
		for ref1, ref2 in combinations(list(self.references.values()), 2):
			if ref1.name == ref2.name:
				log.error(f"{ref1.name} and {ref2.name} have same entry name!")
			if ref1.seq_nt == ref2.seq_nt:
				log.error(f"{ref1.name} and {ref2.name} have same seq!")
			if ref1.segment == 'V' and ref2.segment == 'V' and ref1.is_TRA == ref2.is_TRA:
				try:
					distance = bio.hamming_distance(ref1.seq_nt[-90:], ref2.seq_nt[-90:])
					if distance <= 3:
						log.warn(f"Last 90 bases of {ref1.name} and {ref2.name}: " 
						f"Hamming Distance={distance}", indent=1)
				except Exception: # Sequences with len < 90 will be skipped.
					pass