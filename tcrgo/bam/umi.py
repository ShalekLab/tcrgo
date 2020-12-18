import pysam
import textwrap
from collections import Counter, defaultdict
from itertools import combinations
# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

from typing import List, Dict, Iterator, Set, Tuple, Counter, Union
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile
from .reference import Reference, ReferenceDict
from .read import Read
from .cdr3 import CDR3
import tcrgo.bio as bio
from ..log import Log

log = Log(name=__name__)
log.proceed()

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
					#log.verbose(f"{read[segment]} = {ref_loser}, {top_references}")
					index_top_reference = ref_loser.index(top_references[0])
					#log.verbose(f"{index_top_reference}")
					read[segment] = read[segment][index_top_reference]
					ref_winner = read[segment].reference_name
				else:
					root = read[segment][0].reference_name.split('-')[0]
					if all([candidate.reference_name.split('-')[0] == root \
						for candidate in read[segment][1:]]):
						method = "Root"
						read[segment] = read[segment][0]
						ref_winner = root+'-?'
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
				#if read[segment+"_reference"] is None: 
					#log.verbose(read) # DEBUG

	# TODO: Make disclude_relatives an argument
	def get_top_VJ_reads(self, include_relatives: bool=True) -> List[Read]:
		counts_VJ = Counter()
		reads_VJ = defaultdict(list)
		for read in self.get_reads():
			if include_relatives:
				refs_VJ = (read.top_V_reference.root, read.top_J_reference.root)
			else:
				refs_VJ = (read.top_V_reference.name, read.top_J_reference.name)
			reads_VJ[refs_VJ].append(read)
			counts_VJ[refs_VJ] += 1
		self.top_VJ = max(counts_VJ, key=counts_VJ.get)
		self.count_top_VJ = counts_VJ[self.top_VJ]
		self.frequency_top_VJ = self.count_top_VJ / len(self)
		return reads_VJ[self.top_VJ]

	def count_regions(self):
		for read in self.get_reads(): 
			for subregion in ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC", "UNKN"):
				if read.has_region[subregion]:
					self.counts_region[subregion] += 1

	@log.time
	def select_cdr3(self, reads_VJ: List[Read]):
		log.info(f"UMI {self.sequence}")
		seq_cdr3s = dict()
		seq_counts = Counter()
		counts_cdr3_aa = Counter() #DEBUG
		scores = dict() #DEBUG
		for read in reads_VJ:
			#log.verbose(f"{read.cdr3.seq_nt}, {read.cdr3.seq_aa}, {read.cdr3.score}", indent=1)
			seq_cdr3s[read.cdr3.seq_nt] = read.cdr3
			seq_counts[read.cdr3.seq_nt] += 1
			counts_cdr3_aa[read.cdr3.seq_aa] += 1 #DEBUG
			scores[read.cdr3.seq_aa] = read.cdr3.score #DEBUG
		#log.verbose(seq_counts)
		#log.verbose(counts_cdr3_aa)
		#log.verbose(scores)

		seq_counts_ld = seq_counts.copy()
		for seq1, seq2 in combinations(list(seq_counts_ld.keys()), 2):
			if bio.levenshtein_distance(seq1, seq2) == 1:
				seq_counts_ld[seq1] += 1
				seq_counts_ld[seq2] += 1
		#log.verbose(seq_counts_ld)
		top_count = max(seq_counts_ld.values())
		top_seq_counts = {seq:seq_counts[seq] for seq,count in seq_counts_ld.items() if count == top_count}
		if len(top_seq_counts) == 1:
			self.top_cdr3 = seq_cdr3s[next(iter(top_seq_counts))] # Get CDR3 of the top seq
			self.count_top_cdr3 = top_count
			self.frequency_top_cdr3 = top_count / sum(seq_counts.values())
		else: # TODO: Consider instead retaining ties and doing second pass at barcode level.
			transcript = seq_cdr3s[next(iter(top_seq_counts))].transcript # Only the first transcript is reported.
			self.top_cdr3 = CDR3(bio.get_consensus_sequence(top_seq_counts), None, None, transcript)
			self.top_cdr3.get_translation()
			self.count_top_cdr3 = sum(top_seq_counts.values())
			self.frequency_top_cdr3 = self.count_top_cdr3 / sum(seq_counts.values())
		log.verbose(
			f"WINNER: {self.top_cdr3.seq_nt}, {self.top_cdr3.seq_aa} "
			f"{self.count_top_cdr3}, {self.frequency_top_cdr3}"
		)
			
	#########################################################################
	#	Below is unused or deprecated code
	#########################################################################

	def cdr3_candidates(self, minimum_cdr3s): # self.top_cdr3 saved as a Set() that will be resolved
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

		sum_cdr3_counts = sum(counts_cdr3.values())
		if sum_cdr3_counts >= minimum_cdr3s:
			max_count = 0
			top_cdr3s = None
			for cdr3, count in counts_cdr3.items():
				if len(cdr3) >= 15:
					#if cdr3.startswith(("TGT", "TGC")):
					#if cdr3.endswith("TTT") or cdr3.endswith("TTC"):
					if count > max_count:
						top_cdr3s = {cdr3}
						max_count = count
					elif count == max_count:
						top_cdr3s.add(cdr3)
			if top_cdr3s is not None:
				# TODO: LEVENSHTEIN DISTANCE, any remaining ties keep as a set, resolve on barcode level
				if not self.is_complete_cdr3:
					top_cdr3s = {cdr3+'_' for cdr3 in top_cdr3s}
				self.top_cdr3 = top_cdr3s
				self.count_top_cdr3 = max_count
				self.frequency_top_cdr3 = max_count / sum_cdr3_counts

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
		log.verbose(
			f"For {self.sequence} got {count_complete} complete and "
			f"{count_incomplete} incomplete of {count_total} reads."
		)

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
		log.info(f"{counts_top_V}")
		log.verbose(f"{top_V} {frequency_top_V}")
		log.info(f"{counts_top_J}")
		log.verbose(f"{top_J} {frequency_top_J}")
		log.info(f"{count_total}")
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
		
		log.info(f"{self.counts_top_V}")
		log.verbose(f"{top_V} {frequency_top_V}")
		log.info(f"{self.counts_top_J}")
		log.verbose(f"{top_J} {frequency_top_J}")
		log.info(f"{count_total_V} {count_total_J}")