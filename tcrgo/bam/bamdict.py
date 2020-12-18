# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict
import pysam
import textwrap
import os
import os.path as osp

from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterator, Set, Tuple, DefaultDict
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile
FASTA = pysam.libcfaidx.FastaFile

from .read import Read
from .umi import UMI
from .barcode import Barcode
from .reference import ReferenceDict

from ..log import Log
log = Log(name=__name__)
log.proceed()

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

	#def __del__(self):
	#	self.bam.close()
	#	del self

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
	def build(self, id_queries: List[str], refdict: ReferenceDict):
		log.info("Building BAM Index in memory for fetching VJ queries")
		bam_indexed = pysam.IndexedReads(self.bam)
		bam_indexed.build()
		log.info("Built index in memory for fast retrieval")
		count = 0
		report_interval = len(id_queries) // 20
		for id_query in id_queries:
			barcode, umi, query_name = id_query.split('|') 
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
				log.info(f"Processed {count} reads of {len(id_queries)}")

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

	# TODO: Move to Read file, adjust, and call for each read from reconstruct_tcrs instead
	@log.time
	def reconstruct_cdr3s(self, fasta: FASTA, cdr3_positions: Dict[str, int]):
		for read in self.get_reads():
			read.get_cdr3_positions(cdr3_positions)
			read.ref_seq_V = fasta.fetch(read.top_V.reference_name)
			read.ref_seq_J = fasta.fetch(read.top_J.reference_name)
			read.get_cdr3_sequence()

	@log.time
	def cdr3_statistics(self):
		for umi in self.get_umis():
			umi.cdr3_statistics()

	@log.time
	def write_cdr3_info(self, worker: int, output_path: Path):
		cdr3_info_filename = osp.join(output_path, f"cdr3_info{worker}.tsv")
		if osp.isfile(cdr3_info_filename):
			log.warn(f"Deleting aleady-existing {cdr3_info_filename}")
			os.remove(cdr3_info_filename)
		with open(cdr3_info_filename, 'a') as cdr3_info:
			cdr3_info.write(
				"BC_index\tUMI_index\tBC\tUMI\tnReads\t"
				"topVJ_region\ttopVJ_nReads\ttopVJ_freq\t"
				"CDR3_nt\tCDR3_aa\tCDR3_isProductive\t"
				"CDR3_nReads\tCDR3_freq\tCDR3_stopcodons\t"
				"TCR_nt\tTCR_orf\tTCR_aa\t"
				"TRAV_nReads\tTRAJ_nReads\tTRAC_nReads\t"
				"TRBV_nReads\tTRBJ_nReads\tTRBC_nReads\t"
				"UNKN_nReads\n"
			)
			b = 1
			for barcode_seq, barcode in self.items():
				u = 1
				for umi_seq, umi in barcode.items():
					cdr3_info.write(
						f"{b}\t{u}\t{barcode_seq}\t{umi_seq}\t{len(umi)}\t"
						f"{umi.top_VJ}\t{umi.count_top_VJ}\t{umi.frequency_top_VJ:.3f}\t"
						f"{str(umi.top_cdr3.seq_nt)}\t{str(umi.top_cdr3.seq_aa)}\t{umi.top_cdr3.is_productive}\t"
						f"{umi.count_top_cdr3}\t{umi.frequency_top_cdr3:.3f}\t{len(umi.top_cdr3.transcript.stops)}\t"
						f"{str(umi.top_cdr3.transcript.seq_nt)}\t{umi.top_cdr3.transcript.orf}\t{str(umi.top_cdr3.transcript.seq_aa)}\t"
						f"{umi.counts_region['TRAV']}\t{umi.counts_region['TRAJ']}\t{umi.counts_region['TRAC']}\t"
						f"{umi.counts_region['TRBV']}\t{umi.counts_region['TRBJ']}\t{umi.counts_region['TRBC']}\t"
						f"{umi.counts_region['UNKN']}\n"
					)
					u += 1
				b += 1

	#########################################################################
	#	Below is hacky, deprecated code that will be removed eventually
	#########################################################################

	@log.time
	def write_reads(self, fasta: FASTA, cdr3_positions: Dict[str, int]):
		"""
		This function is not intended for the final pipeline, it's just
		a hack to write all data for some statistic analysis
		"""
		count = 0
		total_reads = len(self.get_reads())
		with open("reads.tsv", 'a') as reads:
			reads.write("QNAME\tBarcode\tUMI\t#Align\tFLAGs\tRegions\tScores\tCIGARs\tMAPQs\tPhreds\tEditDistance\thasVJ\ttopV\ttopJ\tCDR3status\tCDR3seq\n")
			for barcode_seq, barcode in self.items():
				for umi_seq, umi in barcode.items():
					for read_qname, read in umi.items():
						flags = []
						regions = []
						scores = []
						cigars = []
						mapping_quals = []
						phred_scores = []
						edit_distances = []
						has_VJ = False
						top_score_J = 0
						top_score_V = 0
						cur_index = 0
						for alignment in read:
							tags = dict(alignment.tags)
							score = tags["AS"]
							if 'V' in alignment.reference_name:
								if score > top_score_V:
									read.top_V = alignment
									top_score_V = score
									read.top_index_V = cur_index
							elif 'J' in alignment.reference_name:
								if score > top_score_J:
									read.top_J = alignment
									top_score_J = score
									read.top_index_J = cur_index
							flags.append(alignment.flag)
							regions.append(alignment.reference_name)
							scores.append(score)
							cigars.append(alignment.cigarstring)
							mapping_quals.append(alignment.mapping_quality)
							phred_scores.append(tags["UQ"])
							edit_distances.append(tags["NM"])
							cur_index += 1
						if read.top_V is not None and read.top_J is not None:
							has_VJ = True
							read.get_cdr3_positions(cdr3_positions)
							read.ref_seq_V = fasta.fetch(read.top_V.reference_name)
							read.ref_seq_J = fasta.fetch(read.top_J.reference_name)
							read.get_cdr3_sequence()
						reads.write(
							f"{read_qname}\t{barcode_seq}\t{umi_seq}\t{str(len(read))}\t"
							f"{','.join(str(f) for f in flags)}\t"
							f"{','.join(regions)}\t"
							f"{','.join(str(s) for s in scores)}\t"
							f"{','.join(cigars)}\t"
							f"{','.join(str(m) for m in mapping_quals)}\t"
							f"{','.join(str(p) for p in phred_scores)}\t"
							f"{','.join(str(e) for e in edit_distances)}\t"
							f"{has_VJ}\t"
							f"{read.top_index_V}\t"
							f"{read.top_index_J}\t"
							f"{read.cdr3_status}\t{read.cdr3}\n"	
						)
						count += 1 
						if count % 100000 == 0:
							log.info(f"Wrote {count} of {total_reads} unique reads.")
							log.info(f"{read.top_index_V}\t{read.top_index_J}")