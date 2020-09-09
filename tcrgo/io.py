import os.path as osp
import subprocess as sp

from pathlib import Path
from typing import List, Tuple, Dict
import pysam
BAM = pysam.libcalignmentfile.AlignmentFile

from .log import Log
log = Log(name=__name__)
log.proceed()

def sort_and_index(bam: Path) -> Path:
	"""
	Given a path to a BAM, sort the BAM if filename does not end with 'Sorted.bam' or
	'Sort.bam'. Then generate a BAM index (.bai) if it does not exist
	"""
	if not str(bam).endswith(("Sorted.bam", "Sort.bam")):
		bam_basename = osp.splitext(osp.basename(bam))[0]
		bam_sorted_name = f"{bam_basename}Sorted.bam"
		log.info(f"Sorting BAM and creating {bam_sorted_name}")
		pysam.sort('-@', '4', '-o', bam_sorted_name, str(bam))
		bam = bam_sorted_name
	
	bam_index = f"{bam}.bai"
	if osp.isfile(bam_index):
		log.info("BAM index detected.")
	else:
		log.info("No index was detected, building index...")
		pysam.index(bam)
	
	return bam

def determine_data_filetype(data: List[Path]) -> bool:
	"""
	Evaluate file extension and return True for ".bam"
	or False for ".fastq.gz" or ".fastq.gz"
	raise Exception for other cases
	"""
	log.info("Determining the filetype of the sequence data...")
	if len(data) == 1 and str(data[0]).endswith(".bam"):
		return True
	elif len(data) == 2 and all(str(path).endswith((".fastq.gz", ".fastq")) for path in data):
		return False
	else:
		log.error("Please either input one '.bam' path or two '.fastq.gz'/'.fastq' paths.")

# TODO: Redo but with Picard
def fastq_to_bam(reference: Path, data: List[Path]) -> Path:
	"""
	Convert FASTQ to BAM using pysam and return the path to the newly created BAM
	"""
	# Get sample name as the basename without extension for read1
	# sample is used as the name of the created SAM / BAM
	pass
	"""
	sample = osp.splitext(read1)[0]
	output = sp.check_output(
		['bwa', 'mem', fasta, read1, read_tcr]
	).strip().decode()
	sam = f"{sample}TCRalign.sam"
	bam = sam.replace(".bam", "Sort.bam")
	pysam.view(output, '-F', '256', 'h', f'save_stdout={sam}') # Can't use -b flag here?
	# TODO: Only way to convert to BAM? Seems redundant to step in sort_index_bam
	pysam.sort('-o', bam, sam)
	return bam
	"""

# TODO: Should this method return Path or BAM?
def parse_data(data: List[Path]) -> BAM:
	"""
	Read the inputted path(s) to the single-end BAM or paired-end FASTQS
	and determine the file type. Convert FASTQ to BAM if necessary and
	ultimately sort and index the BAM.
	"""
	is_bam = determine_data_filetype(data)
	if is_bam:
		bam = data[0]
	else:
		log.error("FASTQ TO BAM CURRENTLY NOT SUPPORTED") #TODO
		#bam = fastq_to_bam(data)
	bam = sort_and_index(bam)
	log.info(f"Reading {bam}...")
	return pysam.AlignmentFile(bam, "rb")

def output_queries(queries, workers):
	pass

def read_cdr3_file(cdr3_file: str) -> Dict[str, int]:
	cdr3_positions = {}
	with open(cdr3_file, 'r') as cdr3_file:
		for line in cdr3_file:
			line = line.strip().split('\t')
			cdr3_positions[line[0]] = int(line[1])
	return cdr3_positions

"""
	log.info("Reading CDR3bases and FASTA files")
	cdr3_positions = reference.read_cdr3_file(cdr3_file)

	pysam.samtools.faidx(fasta_file)
	fasta = pysam.FastaFile(fasta_file)

	read.get_cdr3_position(read.top_J, cdr3_positions)
	read.get_cdr3_position(read.top_V, cdr3_positions)

	J_seq_ref = fasta.fetch(
		read.top_J.reference_name,
		read.top_J.reference_start,
		read.top_J.reference_end
	)
	V_seq_ref = fasta.fetch(
		read.top_V.reference_name,
		read.top_V.reference_start,
		read.top_V.reference_end
	)
"""
"""
parser.add_argument(
		'-r', "--reference-fasta", 
		#type=argparse.FileType('r', encoding='UTF-8'),
		type=Path,
		required=False,
		help="Path to reference FASTA. Required."
	)
"""