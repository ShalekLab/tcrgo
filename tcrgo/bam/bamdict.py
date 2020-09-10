# TODO: See if can reimplement class using this article:
# https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/ 
#from collections import UserList, UserDict

import pysam
IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile

from typing import List, Dict, Iterator, Set, Tuple, DefaultDict
import textwrap
from collections import Counter, defaultdict

from .read import Read
from .umi import UMI
from .barcode import Barcode

from ..log import Log
log = Log(name=__name__)
log.proceed()

class BAMDict(object):
	"""
	BAMDict[Barcode][UMI][Read].alignments
	"""
	
	def __init__(self, bam: BAM):
		self.barcodes = {}
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

	def __del__(self):
		self.close()
		del self

	def _keys(self) -> List[str]:
		return self.barcodes.keys()

	def _values(self) -> List[Barcode]:
		return self.barcodes.values()

	def items(self) -> List[Tuple[str, Barcode]]:
		return self.barcodes.items()

	def close(self):
		self.bam.close()

	@log.time
	def parse_queries(self) -> Tuple[Set[str], DefaultDict[str, Set[str]]]:
		queries_V = set()
		queries_J = set()
		queries_regions = defaultdict(set)
		for alignment in self.bam.fetch():
			if alignment.is_unmapped:
				log.verbose(f"Read {alignment.query_name} is unmapped, skipping")
				continue
			queries_regions[alignment.query_name].add(alignment.reference_name)
			if 'V' in alignment.reference_name:
				queries_V.add(alignment.query_name)
			elif 'J' in alignment.reference_name:
				queries_J.add(alignment.query_name)
		queries_VJ = queries_V & queries_J
		return queries_VJ, queries_regions

	@log.time		
	def parse_bam(self, query_names: Set[str]): #queries_regions: DefaultDict[str, Set[str]]):
		
		log.info("Building BAM Index for fetching VJ queries")
		bam_indexed = pysam.IndexedReads(self.bam)
		bam_indexed.build()
		log.info("Built")

		count = 0
		for query_name in query_names:
			read = Read(query_name)
			#read.unique_subregions = queries_regions[query_names]
			read.get_alignments(bam_indexed)
			if read.top_V is None or read.top_J is None:
				raise Exception("Couldn't get a top V or J for {query_name}")

			tags = dict(read.top_V.tags)
			barcode = tags["XC"]
			if barcode not in self:
				self[barcode] = Barcode(barcode)
			umi = tags["XU"]
			if umi not in self[barcode]:
				self[barcode][umi] = UMI(umi)
			self[barcode][umi][query_name] = read

			count += 1
			if count % 10000 == 0:
				log.info(f"Processed {count} reads of {len(query_names)}")

	def get_barcode_sequences(self) -> List[str]:
		return self._keys()

	def get_barcodes(self) -> List[Barcode]:
		return self._values()

	def get_umi_sequences(self) -> List[str]:
		umi_seqs = []
		for barcode in self.get_barcodes():
			for umi_seq in barcode.get_umi_sequences():
				umi_seqs.append(umi_seq)

	def get_umis(self) -> List[UMI]:
		umis = []
		for barcode in self.get_barcodes():
			for umi in barcode.get_umis():
				umis.append(umi)
		return umis

	def get_read_query_names(self) -> List[str]:
		read_qnames = []
		umis = self.get_umis()
		for umi in umis:
			for read_qname in umi.get_read_query_names():
				read_qnames.append(read_qname)
		return read_qnames

	def get_reads(self) -> List[Read]:
		reads = []
		umis = self.get_umis()
		for umi in umis:
			for read in umi.get_reads():
				reads.append(read)
		return reads

	def _OLD_parse_bam(self, bam: BAM):
		"""Different approach"""
		for alignment in bam.fetch():
			tags = dict(alignment.tags)
			barcode = tags["XC"]
			umi = tags["XU"]
			read = alignment.query_name
			log.info(f"{barcode} {umi} {read}")

			if alignment.is_unmapped:
				log.verbose(f"Read {alignment} is unmapped, skipping")
				continue
			
			if barcode not in self:
				log.verbose(f"Adding Barcode {barcode}")
				self[barcode] = Barcode(barcode)
			if umi not in self[barcode]:
				log.verbose(f"Adding UMI {umi}")
				self[barcode][umi] = UMI(umi)
			if read not in self[barcode][umi]:
				log.verbose(f"Adding Read {umi}")
				self[barcode][umi][read] = Read(read)
			
			self[barcode][umi][read].append(alignment)
			self[barcode][umi][read].unique_subregions.add(alignment.reference_name)
			log.verbose(f"Current # Alignments: {len(self[barcode][umi][read])}")

	def _discard_single_alignments(self):
		"""Used with the _OLD_parse_bam approach."""
		single_alignment_reads = []
		for barcode_seq, barcode in self.items():
			for umi_seq, umi in barcode.items():
				for read_qname in umi.get_read_query_names():
					if len(self[barcode_seq][umi_seq][read_qname]) == 1:
						single_alignment_reads.append( (barcode_seq, umi_seq, read_qname) )
		for barcode_seq, umi_seq, read_qname in single_alignment_reads:
			del self[barcode_seq][umi_seq][read_qname]



	