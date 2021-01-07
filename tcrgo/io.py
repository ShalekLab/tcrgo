import os
import glob
import random
from collections import Counter, deque, defaultdict
import pandas as pd
import pysam
from .dropseq_tools import fastq_to_bam
from typing import List, Tuple, Dict, Set, DefaultDict, Deque, Iterable
from .log import Log

log = Log("root")
BAM = pysam.libcalignmentfile.AlignmentFile
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
DataFrame = pd.core.frame.DataFrame

def sort_and_index(bam: str, output_path: str) -> str:
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

def parse_data(data: List[str], output_path: str) -> BAM:
	"""
	Read the inputted path(s) to the single-end BAM or paired-end FASTQS
	and determine the file type. Convert FASTQ to BAM if necessary and
	ultimately sort and index the BAM.
	"""
	if is_bam(data):
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
	(Barcode_seq, UMI_seq, QNAME)'
	"""
	queries_V = set()
	queries_J = set()
	queries_all = set()
	count_unmapped = 0
	count_reverse = 0
	for alignment in bam: # .fetch():
		barcode = alignment.get_tag(barcode_tag)
		umi = alignment.get_tag(umi_tag)
		id_query = (barcode, umi, alignment.query_name)
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
	log.info(f"Found and ignored {count_unmapped} unmapped alignments "
		f"and {count_reverse} reverse-mapped alignments.")
	queries_VJ = queries_V & queries_J
	queries_VnotJ = queries_V - queries_J
	queries_JnotV = queries_J - queries_V
	log.info(f"Found and ignored {len(queries_VnotJ)} V-exclusive and "
		f"{len(queries_JnotV)} J-exclusive reads.")
	queries_nonVJ = queries_all - queries_VJ
	return queries_VJ, queries_nonVJ

def filter_queries(queries_VJ: Set[str], num_queries_min: int, \
	num_queries_max: int, seed: int=2020) \
	-> DefaultDict[Tuple[str, str], List[Tuple[str, str, str]]]:
	"""Group queries by BC-UMI combination and apply min and max filters"""
	id_queries = defaultdict(list)
	for id_query in queries_VJ:
		barcode_umi = id_query[:2]
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

def get_partition_queries(w: int, workers: int, max_size: int,
	id_queries: DefaultDict[Tuple[str, str], List[Tuple[str, str, str]]],
	ids_sorted: Deque[Tuple[str, str]],
	id_querycounts: Dict[Tuple[str, str], int]) -> List[str]:
	"""
	Divy up VJ queries as evenly as possible into w files for processing
	by multiple instances of the CDR3 recovery part of the pipeline.
	If I could redo this method I would have it take the deque and a
	DefaultDict[Tuple[str,str], Counter[Tuple[str,str,str], int]] instead.
	"""
	LEFT = 0 # Index of current maximum in ids_sorted
	RIGHT = -1 # Index of current minimum in ids_sorted
	queries = list()
	count = 0
	# TODO: Consider adding: `or # remaining reads < max_size/10`
	if w == workers - 1: # If the last file, just write whatever remains.
		while ids_sorted:
			count += id_querycounts[ids_sorted[LEFT]]
			queries += id_queries[ids_sorted[LEFT]]
			ids_sorted.popleft()
	elif id_querycounts[ids_sorted[LEFT]] >= max_size:
		log.info(f"Sole Barcode-UMI {ids_sorted[LEFT]} of {id_querycounts[ids_sorted[LEFT]]} "
				f"unique reads exceeds max partition size {max_size}.")
		log.warn("Be aware that fewer files may be outputted as a result!")
		count += id_querycounts[ids_sorted[LEFT]]
		queries += id_queries[ids_sorted[LEFT]]
		ids_sorted.popleft()
	else: # current max doesn't below max_size, fill remainder
		while ids_sorted:
			if id_querycounts[ids_sorted[LEFT]] <= max_size - count:
				count += id_querycounts[ids_sorted[LEFT]]
				queries += id_queries[ids_sorted[LEFT]]
				ids_sorted.popleft()
			elif id_querycounts[ids_sorted[RIGHT]] <= max_size - count:
				count += id_querycounts[ids_sorted[RIGHT]]
				queries += id_queries[ids_sorted[RIGHT]]
				ids_sorted.pop()
			else: # No more UMIs can fit into max_size
				break
	return queries

def output_queries(id_queries: Iterable[Tuple[str, str, str]], output_file: str):
	if os.path.isfile(output_file):
		log.warn(f"Deleting pre-existing {output_file}")
		os.remove(output_file)
	with open(output_file, 'w') as query_file:
		for id_query in id_queries:
			query_file.write('\t'.join(id_query) + '\n')

# TODO: consider compressing these files. They can be tens/hundreds of MB
@log.time
def output_grouped_VJ_queries(id_queries: DefaultDict[Tuple[str, str], List[Tuple[str, str, str]]],
	workers: int, output_path: str):
	""" """
	id_querycounts = {k:len(v) for k,v in id_queries.items()} # Convert each value from list to len(value)
	ids_sorted = deque([k for k,v in sorted(id_querycounts.items(), key=lambda item: item[1], reverse=True)]) # Sort keys in descending order
	
	num_queries = sum(id_querycounts.values())
	partition_size = num_queries // workers
	remainder = num_queries % partition_size
	max_size = partition_size + remainder
	log.info(f"Distributing {num_queries} queries belonging to {len(ids_sorted)} "
			f"Barcode-UMI pairs somewhat evenly across {workers} file(s)...")
	log.sep('-', width=50)
	for w in range(workers):	
		if not ids_sorted:
			log.warn(f"No reads remain, not writing further files.")
			break
		queries = get_partition_queries(w, workers, max_size, id_queries, ids_sorted, id_querycounts)
		filename = f"queries{w+1}.tsv"
		output_queries(queries, os.path.join(output_path, filename))
		log.info(f"Wrote {len(queries)} unique reads to {filename}.")
		log.sep("-", width=50)
	if ids_sorted:
		log.error(f"Ended with {len(ids_sorted)} still left!")

def read_queries(queries_filename: str) -> List[Tuple[str, str, str]]:
	log.info(f"Reading {queries_filename}.")
	id_queries = open(queries_filename, 'r').read().splitlines()
	return [id_query.split('\t') for id_query in id_queries]

def get_worker_id(query_list: str) -> int:
	return int(
		os.path.basename(query_list) \
			.replace("queries", '') \
			.replace(".tsv", '')
	)

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

def read_cdr3_info(workers: Iterable[int], input_path: str,
	string_index: bool=False) -> DataFrame:
	"""
	Read the CDR3 info files for the range of workers, concat into dataframe,
	then sort alphabetically by BC and UMI.
	By default, will add numerical indices to the dataframe for BC and UMI.
	If string_index: no numerical indices added but BC repeats are replaced
	with blank cells to make reading by BC easier.
	"""
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

def read_tiebreaks_alignments(workers, input_path: str) -> DataFrame:
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

def write_dataframe(df: DataFrame, output_path: str, filename: str):
	filename_path = os.path.join(output_path, filename)
	log.info(f"Writing aggregated dataframe to file {filename_path}.")
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	if os.path.isfile(filename_path):
		log.warn(f"Deleting already existing {filename_path}.")
		os.remove(filename_path)
	df.to_csv(filename_path, sep='\t', header=True, index=False)
	log.info(f"Wrote {filename_path}")