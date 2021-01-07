import pysam
import textwrap
import os
import os.path as osp
from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterator, Set, Tuple, DefaultDict
from .read import Read
from .umi import UMI
from .barcode import Barcode
from .reference import ReferenceDict
from ..log import Log

log = Log("root")
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile
FASTA = pysam.libcfaidx.FastaFile

class BAMDict(object):
	"""
	BAMDict[Barcode][UMI][Read].alignments
	"""
	
	def __init__(self, bam: BAM):
		self.barcodes = dict()
		self.bam = bam

	def __len__(self) -> int:
		return len(self.barcodes)

	def __getitem__(self, key: str) -> Barcode:
		return self.barcodes[key]

	def __setitem__(self, key, value):
		self.barcodes[key] = value

	def __delitem__(self, key):
		del self.barcodes[key]

	def __iter__(self) -> Iterator:
		return iter(self.barcodes)

	def __contains__(self, item) -> bool:
		return item in self.barcodes

	def __missing__(self, key):
		raise KeyError("No Barcode in BAMDict")

	def _keys(self) -> List[str]:
		return self.barcodes.keys()

	def _values(self) -> List[Barcode]:
		return self.barcodes.values()

	def items(self) -> List[Tuple[str, Barcode]]:
		return self.barcodes.items()

	def close(self):
		self.bam.close()

	def get_barcode_sequences(self) -> List[str]:
		return self._keys()

	def get_barcodes(self) -> List[Barcode]:
		return self._values()

	def get_umi_sequences(self) -> List[str]:
		umi_seqs = []
		for barcode in self.get_barcodes():
			umi_seqs += barcode.get_umi_sequences()

	def get_umis(self) -> List[UMI]:
		umis = []
		for barcode in self.get_barcodes():
			umis += barcode.get_umis()
		return umis

	def get_read_query_names(self) -> List[str]:
		read_qnames = []
		for umi in self.get_umis():
			read_qnames += self.get_read_query_names()
		return read_qnames

	def get_reads(self) -> List[Read]:
		reads = [] 
		for umi in self.get_umis():
			reads += umi.get_reads()
		return reads

	@log.time
	def build(self, id_queries: List[Tuple[str, str, str]], refdict: ReferenceDict):
		log.info("Building BAM Index in memory for fetching VJ queries")
		bam_indexed = pysam.IndexedReads(self.bam)
		bam_indexed.build()
		log.info("Built index in memory for fast retrieval")
		count = 0
		report_interval = max(1, len(id_queries) // 25)
		for id_query in id_queries:
			barcode, umi, query_name = id_query
			read = Read(query_name)
			read.parse_alignments(bam_indexed, refdict)
			if read.top_V is None or read.top_J is None:
				log.error(f"Should have gotten a top V and J for {query_name} but did not!")
			if barcode not in self:
				self[barcode] = Barcode(barcode)
			if umi not in self[barcode]:
				self[barcode][umi] = UMI(umi)
			self[barcode][umi][query_name] = read
			count += 1
			if count % report_interval == 0:
				log.info(f"Stored top alignments for {count} reads ({(count/len(id_queries)):.0%}).")
		del bam_indexed

	def write_tiebreaks_alignments(self, w, output_path):
		ties_aggregated = Counter()
		for read in self.get_reads():
			ties_aggregated += read.ties
		filename = osp.join(output_path, f"tiebreaks_alignments{w}.tsv")
		with open(filename, 'w') as tiebreaks:
			tiebreaks.write("Winner\tLoser\tMethod\tCount\n")
			for case, count in ties_aggregated.items():
				winner, loser, method = case
				tiebreaks.write(f"{winner}\t{loser}\t{method}\t{count}\n")		

	'''
	# TODO: Move to Read file, adjust, and call for each read from reconstruct_tcrs instead
	@log.time
	def reconstruct_cdr3s(self, fasta: FASTA, cdr3_positions: Dict[str, int]):
		for read in self.get_reads():
			read.get_cdr3_positions(cdr3_positions)
			read.ref_seq_V = fasta.fetch(read.top_V.reference_name)
			read.ref_seq_J = fasta.fetch(read.top_J.reference_name)
			read.get_cdr3_sequence()
	'''

	@log.time
	def write_cdr3_info(self, worker: int, output_path: Path):
		cdr3_info_filename = osp.join(output_path, f"cdr3_info{worker}.tsv")
		if osp.isfile(cdr3_info_filename):
			log.warn(f"Deleting aleady-existing {cdr3_info_filename}")
			os.remove(cdr3_info_filename)
		with open(cdr3_info_filename, 'a') as cdr3_info:
			cdr3_info.write('\t'.join([
				"BC_index", "UMI_index", "BC", "UMI", "nReads",
				"topVJ_region", "topVJ_nReads", "topVJ_freq",
				"CDR3_nt", "CDR3_aa", "CDR3_isProductive",
				"CDR3_nReads", "CDR3_freq", "CDR3_stopcodons",
				"TCR_nt", "TCR_orf", "TCR_aa",
				"TRAV_nReads", "TRAJ_nReads", "TRAC_nReads",
				"TRBV_nReads", "TRBJ_nReads", "TRBC_nReads",
				"UNKN_nReads\n"
			]))
			b = 1
			for barcode_seq, barcode in self.items():
				u = 1
				for umi_seq, umi in barcode.items():
					cdr3_info.write('\t'.join([str(elem) for elem in \
						[b, u, barcode_seq, umi_seq, len(umi),
						umi.top_VJ, umi.count_top_VJ, f"{umi.frequency_top_VJ:.3f}",
						umi.top_cdr3.seq_nt, umi.top_cdr3.seq_aa, umi.top_cdr3.is_productive,
						umi.count_top_cdr3, f"{umi.frequency_top_cdr3:.3f}", len(umi.top_cdr3.transcript.stops),
						umi.top_cdr3.transcript.seq_nt, umi.top_cdr3.transcript.orf, umi.top_cdr3.transcript.seq_aa,
						umi.counts_region['TRAV'], umi.counts_region['TRAJ'], umi.counts_region['TRAC'],
						umi.counts_region['TRBV'], umi.counts_region['TRBJ'], umi.counts_region['TRBC'],
						f"{umi.counts_region['UNKN']}\n"]
					]))
					u += 1
				b += 1
		