
import os
import os.path as osp
import subprocess as sp

from pathlib import Path
from typing import List, Tuple, Dict, Set
import pysam
BAM = pysam.libcalignmentfile.AlignmentFile

from .log import Log
log = Log(name=__name__)
log.proceed()

def sort_and_index(bam: Path, output_path: str) -> Path:
	"""
	Given a path to a BAM, sort the BAM if filename does not end with '_sorted.bam'.
	Output to sorted BAM to output_path and generate a BAM index (.bai) if it does not exist
	If a sorted BAM already exists at output_path, use that instead
	"""
	if not osp.exists(output_path):
		os.mkdir(output_path)
	bam_basename = osp.splitext(osp.basename(bam))[0]
	bam_sorted_name = f"{bam_basename}_sorted.bam"
	bam_sorted_path = osp.join(output_path, bam_sorted_name)

	if str(bam).endswith("_sorted.bam"):
		bam_sorted_path = str(bam)
	elif osp.isfile(bam_sorted_path):
		log.warn("Sorted BAM already found at the output path! Using this BAM instead!")
	else:
		log.info(f"Sorting BAM by coordinates, outputting as {bam_sorted_name}")
		pysam.sort('-@', '4', '-o', bam_sorted_path, str(bam))

	bam_index = f"{bam_sorted_path}.bai"
	if osp.isfile(bam_index):
		log.info("BAM index detected.")
	else:
		log.info("No index was detected, building index...")
		pysam.index(bam_sorted_path)
	log.info(f"Reading {bam_sorted_path}...")
	return pysam.AlignmentFile(bam_sorted_path, "rb")

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
def parse_data(data: List[Path], output_path: Path) -> BAM:
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
	bam = sort_and_index(bam, output_path)
	return bam

# TODO: consider compressing these files. They can be tens/hundreds of MB
def output_queries(queries: Set[str], output_path: Path, workers: int):
	queries = list(queries)
	p = len(queries) // workers # partition size
	r = len(queries) % p # remainder
	for i in range(workers):
		queries_filename = osp.join(output_path, f"queries{i+1}.txt")
		if osp.isfile(queries_filename):
			log.warn(f"Deleting aleady-existing {queries_filename}")
			os.remove(queries_filename)
		with open(queries_filename, 'w') as query_list:
			query_list.write('\n'.join(queries[i*p + min(i, r) : (i+1)*p + min(i+1, r)]))

def input_queries(output_path: Path, worker: int) -> Set[str]:
	queries_filename = osp.join(output_path, f"queries{worker}.txt")
	return set(open(queries_filename, 'r').read().splitlines())

def read_cdr3_file(cdr3_file: str) -> Dict[str, int]:
	cdr3_positions = {}
	with open(cdr3_file, 'r') as cdr3_file:
		for line in cdr3_file:
			line = line.strip().split('\t')
			cdr3_positions[line[0]] = int(line[1])
	return cdr3_positions

"""
parser.add_argument(
		'-r', "--reference-fasta", 
		#type=argparse.FileType('r', encoding='UTF-8'),
		type=Path,
		required=False,
		help="Path to reference FASTA. Required."
	)
"""