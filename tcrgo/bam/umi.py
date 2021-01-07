import pysam
import textwrap
from collections import Counter, defaultdict
from itertools import combinations
from typing import List, Dict, Iterator, Set, Tuple, Counter, Union
from .reference import Reference, ReferenceDict
from .read import Read
from .cdr3 import CDR3
import tcrgo.bio as bio
from ..log import Log

log = Log("root")
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

class UMI(object):
	def __init__(self, sequence: str):
		self.sequence = sequence
		self.reads = dict()
		self.counts_region = Counter()

		self.top_VJ = None
		self.count_top_VJ = 0
		self.frequency_top_VJ = 0
		
		self.reads_top_VJ = None # TODO: Is necessary to store in UMI?
		self.top_cdr3 = None
		self.count_top_cdr3 = 0
		self.frequency_top_cdr3 = 0

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

	@log.time
	def resolve_alignment_identities(self, refdict: ReferenceDict):
		identity_V = Counter()
		identity_J = Counter()
		unresolved_V = list()
		unresolved_J = list()
		for read in self.get_reads():
			for segment, identity, unresolved in ( \
				("top_V", identity_V, unresolved_V), \
				("top_J", identity_J, unresolved_J)):
				if isinstance(read[segment], list):
					unresolved.append(read)
					for candidate in read[segment]:
						identity[candidate.reference_name] += 1
				else:
					identity[read[segment].reference_name] += 1
					read[segment+"_reference"] = refdict[read[segment].reference_name]
		for segment, identity, unresolved in ( \
				("top_V", identity_V, unresolved_V), \
				("top_J", identity_J, unresolved_J)):
			for read in unresolved:
				ref_loser = tuple([alignment.reference_name for alignment in read[segment]])
				candidates = {ref:identity[ref] for ref in ref_loser}
				top_reference_count = max(candidates.values())
				top_references = [ref for ref,c in candidates.items() if c == top_reference_count]
				if len(top_references) == 1:
					method = "UMI"
					index_top_reference = ref_loser.index(top_references[0])
					read[segment] = read[segment][index_top_reference]
					ref_winner = read[segment].reference_name
				else:
					root = read[segment][0].reference_name.split('-')[0]
					if all([candidate.reference_name.split('-')[0] == root \
						for candidate in read[segment][1:]]):
						method = "Root"
						read[segment] = read[segment][0]
						ref_winner = root + '-?'
					else:
						method = "Alphabetical"
						read[segment] = read[segment][0]
						ref_winner = read[segment].reference_name
				if ref_winner in refdict:
					read[segment+"_reference"] = refdict[ref_winner]
				else:
					ref_closest = refdict[read[segment].reference_name]
					read[segment+"_reference"] = Reference(
						ref_winner, ref_closest.seq_nt, 
						ref_closest.cdr3_positions, ref_closest.orf
					)
				read.ties[(ref_winner, ref_loser, method)] += 1

	def get_top_VJ_reads(self, exclude_relatives: bool=False) -> List[Read]:
		"""
		Return the Reads that correspond to the highest counted VJ combination.
		If not exclude_relatives (default), reads with references matching the
		top V or J roots (TRAV11 versus TRAV11-1, -2, etc.) are included
		If exclude_relatives, reads must match the specific variant combination
		exactly to be included! If top variant combo is (TRAV11-1, TRAJ33), 
		only reads with top V and J that match these two references will be 
		included for CDR3 recovery.
		Regardless, the specific variant combination is reported for this UMI.
		"""
		counts_VJ = Counter()
		counts_VJ_root = Counter()
		reads_VJ = defaultdict(list)
		for read in self.get_reads():
			refs_VJ = (read.top_V_reference.name, read.top_J_reference.name)
			counts_VJ[refs_VJ] += 1
			if not exclude_relatives:
				refs_VJ = (read.top_V_reference.root, read.top_J_reference.root)
				counts_VJ_root[refs_VJ] += 1
			reads_VJ[refs_VJ].append(read)
		self.top_VJ = max(counts_VJ, key=counts_VJ.get) # Always report top variant
		if not exclude_relatives:
			top_VJ_root = max(counts_VJ_root, key=counts_VJ_root.get)
			self.count_top_VJ = counts_VJ_root[top_VJ_root]
			reads = reads_VJ[top_VJ_root]
		else:
			self.count_top_VJ = counts_VJ[self.top_VJ]
			reads = reads_VJ[self.top_VJ]
		self.frequency_top_VJ = self.count_top_VJ / len(self)
		return reads

	def count_regions(self):
		for read in self.get_reads(): 
			for subregion in ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC", "UNKN"):
				if read.has_region[subregion]:
					self.counts_region[subregion] += 1

	@log.time
	def select_cdr3(self, reads_VJ: List[Read]):
		log.verbose(f"UMI {self.sequence}")
		seq_cdr3s = dict()
		seq_counts = Counter()
		for read in reads_VJ:
			seq_cdr3s[read.cdr3.seq_nt] = read.cdr3
			seq_counts[read.cdr3.seq_nt] += 1

		seq_counts_ld = seq_counts.copy()
		if len(seq_counts_ld.keys()) <= 5:
			distance_func = bio.levenshtein_distance
		else:
			distance_func = bio.hamming_distance
		for seq1, seq2 in combinations(list(seq_counts_ld.keys()), 2):
			try:
				if distance_func(seq1, seq2) == 1:
					seq_counts_ld[seq1] += 1
					seq_counts_ld[seq2] += 1
			except ValueError:
				pass # If hamming distance and unequal lengths, just move on.
		top_count = max(seq_counts_ld.values())
		top_seq_counts = {
			seq:seq_counts[seq] for seq,count \
				in seq_counts_ld.items() if count == top_count
		}
		if len(top_seq_counts) == 1:
			self.top_cdr3 = seq_cdr3s[next(iter(top_seq_counts))] # Get CDR3 of the top seq
			self.count_top_cdr3 = top_count
			self.frequency_top_cdr3 = top_count / sum(seq_counts.values())
		else: # TODO: Consider instead retaining ties and doing second pass at barcode level.
			seq_nt = bio.get_consensus_sequence(top_seq_counts)
			self.top_cdr3 = CDR3(
				seq_nt=seq_nt, 
				start=0,
				end=len(seq_nt) - 1, 
				transcript=seq_cdr3s[next(iter(top_seq_counts))].transcript
			) # First transcript is chosen to represent
			self.top_cdr3.get_translation()
			self.count_top_cdr3 = sum(top_seq_counts.values())
			self.frequency_top_cdr3 = self.count_top_cdr3 / sum(seq_counts.values())
		log.verbose(
			f"WINNER: {self.top_cdr3.seq_nt}, {self.top_cdr3.seq_aa} "
			f"{self.count_top_cdr3}, {self.frequency_top_cdr3}"
		)
