import os
import glob
import random
from collections import Counter, deque, defaultdict
import pandas as pd
import pysam
from .dropseq_tools import fastq_to_bam
from .log import Log

BAM = pysam.libcalignmentfile.AlignmentFile
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
DataFrame = pd.core.frame.DataFrame
from pathlib import Path
from typing import List, Tuple, Dict, Set, DefaultDict, Deque

log = Log(name=__name__)
log.proceed()

def sort_and_index(bam: Path, output_path: str) -> Path:
	"""
	Given a path to a BAM, sort the BAM if filename does not end with 'sorted.bam'.
	Output to sorted BAM to output_path and generate a BAM index (.bai) if it does not exist
	If a sorted BAM already exists at output_path, use that instead
	"""
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	bam_basename = os.path.splitext(os.path.basename(bam))[0]
	bam_sorted_name = f"{bam_basename}sorted.bam"
	bam_sorted_path = os.path.join(output_path, bam_sorted_name)

	if str(bam).endswith("sorted.bam"):
		bam_sorted_path = str(bam)
	elif os.path.isfile(bam_sorted_path):
		log.warn("Sorted BAM already found at the output path! Using this BAM instead!")
	else:
		log.info(f"Sorting BAM by coordinates, outputting as {bam_sorted_name}")
		pysam.sort('-@', '4', '-o', bam_sorted_path, str(bam))

	bam_index = f"{bam_sorted_path}.bai"
	if os.path.isfile(bam_index):
		log.info("BAM index detected.")
	else:
		log.info("No index was detected, building index...")
		pysam.index(bam_sorted_path)
	log.info(f"Reading {bam_sorted_path}...")
	return pysam.AlignmentFile(bam_sorted_path, "rb")

def is_bam(data: List[str]) -> bool:
	"""
	Evaluate file extension and return True for ".bam"
	or False for ".fastq.gz" or ".fastq.gz"
	raise Exception for other cases
	"""
	log.info("Determining the filetype of the sequence data...")
	if len(data) == 1 and data[0].endswith(".bam"):
		return True
	elif all(path.endswith((".fastq.gz", ".fastq")) for path in data):
		return False
	else:
		log.error("Please either input one '.bam' path or two '.fastq.gz'/'.fastq' paths.")

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
		log.error("Please run the alignment script to convert paired end FASTQs to an aligned BAM.")
	bam = sort_and_index(bam, output_path)
	return bam

@log.time
def parse_queries(bam: BAM, barcode_tag: str="CR", umi_tag: str="RX") -> Tuple[Set[str], Set[str]]:
	"""
	Parses BAM line by line and finds intersection of queries containing V 
	and queries containing J. Stores these queries as
	'Barcode_seq|UMI_seq|QNAME'
	"""
	queries_V = set()
	queries_J = set()
	queries_all = set()
	count_unmapped = 0
	count_reverse = 0
	for alignment in bam: # .fetch():
		barcode = alignment.get_tag(barcode_tag)
		umi = alignment.get_tag(umi_tag)
		id_query = f"{barcode}|{umi}|{alignment.query_name}"
		queries_all.add(id_query)
		if alignment.is_unmapped:
			count_unmapped += 1
			continue
		if alignment.is_reverse:
			count_reverse += 1
			continue
		if "TRAV" in alignment.reference_name or "TRBV" in alignment.reference_name:
			if alignment.reference_end >= bam.get_reference_length(alignment.reference_name) - 20:
				queries_V.add(id_query)
		elif "TRAJ" in alignment.reference_name or "TRBJ" in alignment.reference_name:
			queries_J.add(id_query)
	log.info(f"Found and ignored {count_unmapped} unmapped reads "
		f"and {count_reverse} reverse mapped reads.")
	queries_VJ = queries_V & queries_J
	queries_VnotJ = queries_V - queries_J
	queries_JnotV = queries_J - queries_V
	log.info(f"Found and ignored {len(queries_VnotJ)} V-exclusive and "
		f"{len(queries_JnotV)} J-exclusive reads.")
	queries_nonVJ = queries_all - queries_VJ
	return queries_VJ, queries_nonVJ

def filter_queries(queries_VJ: Set[str], num_queries_min: int, \
	num_queries_max: int, seed: int=2020) -> DefaultDict[str, List[str]]:
	""" """
	id_queries = defaultdict(list)
	for id_query in queries_VJ:
		barcode_umi = '|'.join(id_query.split('|')[:2])
		id_queries[barcode_umi].append(id_query)
	ids_to_drop = list()
	random.seed(seed)
	for barcode_umi, queries in id_queries.items():
		num_queries = len(queries)
		if num_queries < num_queries_min:
			ids_to_drop.append(barcode_umi)
		elif num_queries > num_queries_max:
			id_queries[barcode_umi] = random.sample(queries, num_queries_max)
	log.info(f"Dropping {len(ids_to_drop)} queries below the minimum threshold")
	for barcode_umi in ids_to_drop:
		del id_queries[barcode_umi]
	return id_queries

def get_partition_queries(w: int, workers:int, max_size: int, \
	id_queries: DefaultDict[str, List[str]], ids_sorted: Deque[str], \
	id_querycounts: Dict[str, int]) -> List[str]:
	"""
	Divy up VJ queries as evenly as possible into w files for processing
	by multiple instances of the CDR3 recovery part of the pipeline.
	"""
	LEFT = 0 # Index of current maximum in ids_sorted
	RIGHT = -1 # Index of current minimum in ids_sorted
	to_write = list()
	count = 0
	# TODO: Consider adding: `or # remaining reads < max_size/10`
	if w == workers - 1: # If the last file, just write whatever remains.
		while ids_sorted:
			count += id_querycounts[ids_sorted[LEFT]]
			to_write += id_queries[ids_sorted[LEFT]]
			ids_sorted.popleft()
	elif id_querycounts[ids_sorted[LEFT]] >= max_size:
		log.info(f"Sole Barcode|UMI {ids_sorted[LEFT]} of {id_querycounts[ids_sorted[LEFT]]} "
				f"unique reads exceeds max partition size {max_size}.")
		log.warn("Be aware that fewer files may be outputted as a result!")
		count += id_querycounts[ids_sorted[LEFT]]
		to_write += id_queries[ids_sorted[LEFT]]
		ids_sorted.popleft()
	else: # current max doesn't below max_size, fill remainder
		while ids_sorted:
			if id_querycounts[ids_sorted[LEFT]] <= max_size - count:
				count += id_querycounts[ids_sorted[LEFT]]
				to_write += id_queries[ids_sorted[LEFT]]
				ids_sorted.popleft()
			elif id_querycounts[ids_sorted[RIGHT]] <= max_size - count:
				count += id_querycounts[ids_sorted[RIGHT]]
				to_write += id_queries[ids_sorted[RIGHT]]
				ids_sorted.pop()
			else: # No more UMIs can fit into max_size
				break
	return to_write

# TODO: consider compressing these files. They can be tens/hundreds of MB
@log.time
def output_grouped_VJ_queries(id_queries: DefaultDict[str, List[str]], workers: int, output_path: Path):
	id_querycounts = {k:len(v) for k,v in id_queries.items()} # Convert each value from list to len(list)
	ids_sorted = deque([k for k,v in sorted(id_querycounts.items(), key=lambda item: item[1], reverse=True)]) # Sort keys in descending order
	
	num_queries = sum(id_querycounts.values())
	partition_size = num_queries // workers
	remainder = num_queries % partition_size
	max_size = partition_size + remainder
	log.info(f"Distributing {num_queries} queries belonging to {len(ids_sorted)} "
			f"Barcode-UMI pairs somewhat evenly across {workers} file(s).")
	log.sep('-', width=50)

	for w in range(workers):	
		queries_filename = os.path.join(output_path, f"queries{w+1}.txt")
		if not ids_sorted:
			log.warn(f"No reads remain, not writing file '{queries_filename}'!")
			continue
		if os.path.isfile(queries_filename):
			log.warn(f"Deleting aleady-existing {queries_filename}")
			os.remove(queries_filename)
		to_write = get_partition_queries(w, workers, max_size, id_queries, ids_sorted, id_querycounts)
		with open(queries_filename, 'w') as query_list:
			query_list.write('\n'.join(to_write)+'\n')
		log.info(f"Wrote {len(to_write)} unique reads to file {w+1}")
		log.sep("-", width=50)
	if ids_sorted:
		log.error(f"Ended with {len(ids_sorted)} still left!")

def output_nonVJ(queries_nonVJ: Set[str], output_path: Path):
	filename = os.path.join(output_path, f"queries_nonVJ.txt")
	with open(filename, 'w') as nonVJ_file:
		nonVJ_file.write('\n'.join(queries_nonVJ)+'\n')

def read_id_queries(output_path: Path, worker: int) -> List[str]:
	queries_filename = os.path.join(output_path, f"queries{worker}.txt")
	log.info(f"Reading {queries_filename}.")
	return open(queries_filename, 'r').read().splitlines()

def read_query_list(query_list: str) -> List[str]:
	log.info(f"Reading {query_list}.")
	return open(query_list, 'r').read().splitlines()

def get_worker_id(query_list: str) -> int:
	return int(os.path.basename(query_list).replace("queries", '').replace(".txt", ''))

def parse_fasta(fasta: str) -> List[Tuple[str, str]]:
	entries = list()
	with open(fasta, 'r') as fasta:
		for line in fasta:
			if line.startswith('>'):
				entry = line.replace('>', '').strip()
				seq = fasta.readline().strip()
				entries.append((entry, seq))
	if not entries:
		log.error("Found no entries beginning with '>' in the FASTA")
	return entries

def read_cdr3_positions(cdr3_file: str) -> Dict[str, int]:
	cdr3_positions = {}
	with open(cdr3_file, 'r') as cdr3_file:
		for line in cdr3_file:
			line = line.strip().split('\t')
			cdr3_positions[line[0]] = int(line[1])
	return cdr3_positions

def list_cdr3_files(input_directory: str) -> List[int]:
	cdr3_infos = glob.glob(os.path.join(input_directory, "cdr3_info*.tsv"))
	worker_ids = list()
	for cdr3_info in cdr3_infos:
		w = os.path.basename(cdr3_info).replace("cdr3_info", '').replace(".tsv", '')
		worker_ids.append(int(w))
	return sorted(worker_ids)

def read_cdr3_info(workers, input_path: Path, string_index: bool=False) -> DataFrame:
	for w in workers:
		cdr3_info_filename = os.path.join(input_path, f"cdr3_info{w}.tsv")
		w_cdr3_info = pd.read_csv(cdr3_info_filename, sep='\t', header=0, index_col=["BC", "UMI"])
		if w == 1:
			df = w_cdr3_info
		else:
			df = pd.concat([df, w_cdr3_info])
	
	log.info("Sorting aggregated DataFrame by 'BC' and 'UMI'.")
	df = df.sort_values(["BC", "UMI"])
	df.drop(["BC_index", "UMI_index"], axis=1, inplace=True)
	if string_index:
		mask = df.index.get_level_values("BC").to_series().duplicated()
		df = df.reset_index()
		df.loc[mask.values,['BC']] = ''
	else:
		umi_index = df.groupby(level="BC").cumcount()+1
		umi_index = umi_index.reset_index(drop=True)
		counts = pd.DataFrame(df.groupby("BC").size()).rename(columns={0:"count"}).reset_index(drop=True)
		counts.index += 1
		bc_index = counts.index.to_series().repeat(counts["count"]).reset_index(drop=True)
		df = df.reset_index()
		df.insert(0, "UMI_index", umi_index)
		df.insert(0, "BC_index", bc_index)	
	return df

def read_tiebreaks_alignments(workers, input_path: Path) -> DataFrame:
	for w in workers:
		tiebreaks_alignments = os.path.join(input_path, f"tiebreaks_alignments{w}.tsv")
		w_tiebreaks_alignments = pd.read_csv(tiebreaks_alignments, sep='\t', header=0)
		if w == 1:
			df = w_tiebreaks_alignments
		else:
			df = df.merge(w_tiebreaks_alignments, how="outer", on=["Winner", "Loser", "Method"])
			df = df.fillna(0)
			df["Count"] = (df["Count_x"] + df["Count_y"]).astype(int)
			df = df.drop(["Count_x", "Count_y"], axis=1)
	df = df.sort_values(["Winner", "Loser", "Method"])
	return df

def write_dataframe(df: DataFrame, output_path: Path, filename: str):
	filename_path = os.path.join(output_path, filename)
	log.info(f"Writing aggregated dataframe to file {filename_path}.")
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	if os.path.isfile(filename_path):
		log.warn(f"Deleting already existing {filename_path}.")
		os.remove(filename_path)
	df.to_csv(filename_path, sep='\t', header=True, index=False)
	log.info(f"Wrote {filename_path}")

###################################################################################
#	Deprecated
###################################################################################

def output_queries(queries: Set[str], output_path: Path, workers: int):
	queries = list(queries)
	p = len(queries) // workers # partition size
	r = len(queries) % p # remainder
	for i in range(workers):
		queries_filename = os.path.join(output_path, f"queries{i+1}.txt")
		if os.path.isfile(queries_filename):
			log.warn(f"Deleting aleady-existing {queries_filename}")
			os.remove(queries_filename)
		with open(queries_filename, 'w') as query_list:
			query_list.write('\n'.join(queries[i*p + min(i, r) : (i+1)*p + min(i+1, r)]))