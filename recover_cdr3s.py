"""
This is the second script of the python pipeline. 
It is meant to take the filtered reads, find the top V and J alignments,
reconstruct the CDR3 sequences and select the most representative CDR3.
Run using python -m recover_cdr3s <args>

Requirements:
	Python >3.8.5, samtools, pysam, biopython, pandas
"""
import argparse
import os.path
from textwrap import dedent, indent
import pysam

from typing import List, Dict, Iterator, Set
from pathlib import Path

from tcrgo.bam import BAMDict, ReferenceDict
import tcrgo.io as io
from tcrgo import Log
log = Log("root")

def main(args):
	log.init(args.verbosity)
	log.info("Parsing input data...")
	bam = io.sort_and_index(args.bam, args.output_path)
	
	log.info("Reading CDR3 positions file and FASTA file...")
	cdr3_positions = None
	if args.cdr3_positions_file is not None:
		cdr3_positions = io.read_cdr3_positions(args.cdr3_positions_file)
	log.verbose(f"{cdr3_positions}")
	entries = io.parse_fasta(args.fasta)
	refdict = ReferenceDict()
	refdict.build(entries, cdr3_positions, args.zero_indexed)
	refdict.verify()

	log.info("Fetching filtered queries containing alignments to both V and J segments.")
	if args.query_list is not None:
		queries_filename = args.query_list
		worker = io.get_worker_id(args.query_list)
	else:
		worker = args.worker
		queries_filename = os.path.join(args.input_path, f"queries{worker}.tsv")
	id_queries = io.read_queries(queries_filename)

	bamdict = BAMDict(bam)
	log.info("Finding the top V and J alignments for each query")
	# TODO: Consider splitting fully-mapped V/J reads and partially-mapped V/J reads into different camps
	# Perform de Bruijn graph sequence assembly on partially mapped reads? At least see what they're like.
	bamdict.build(id_queries, refdict)
	log.info("Finished retrieving top V and J alignments for each query.")
	num_umis = len(bamdict.get_umis())
	log.info(f"Got {len(bamdict.get_barcodes())} barcode(s), "
			f"{num_umis} UMI(s) and "
			f"{len(bamdict.get_reads())} read(s) total")
	
	log.info("For each UMI, resolving V and J alignment identities "
	 	"and selecting the most representative CDR3 sequences...")
	count = 0
	count_filtered = 0
	report_interval = max(1, num_umis // 50)
	for barcode in bamdict.get_barcodes():
		for umi in barcode.get_umis():
			umi.count_regions()
			umi.resolve_alignment_identities(refdict)
			reads_VJ = umi.get_top_VJ_reads(args.exclude_relatives)
			if len(reads_VJ) >= args.minimum_cdr3s \
				and umi.frequency_top_VJ >= args.minimum_frequency:
				for read in reads_VJ:
					read.get_cdr3()
				umi.select_cdr3(reads_VJ)
			else:
				count_filtered += 1
			count += 1
			if count % report_interval == 0:
				log.info(f"Processed {count} UMIs ({(count / num_umis):.0%}).")

	log.info("Finished recovering CDR3 sequences for all VJ reads "
		 "and recovered top CDR3 for each UMI.")
	log.info(f"{count_filtered} BCUMIs were filtered out due to top VJ "
		"read number and frequency thresholds set by the program's arguments.")
	log.info(f"Writing the remaining {(count - count_filtered)} BCUMIs to file.")
	bamdict.write_cdr3_info(worker, args.output_path)
	log.info("Writing reference info for alignment tiebreaks.")
	bamdict.write_tiebreaks_alignments(worker, args.output_path)
	bamdict.close()
	log.success("DONE")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		prog="python -m recover_cdr3s", # TODO: Would be good to make a command-line alias and set that here
		usage="",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=dedent(
			f"""
			+{'='*78}+
			|{'{:^78}'.format('SEQ-WELL TCR RECOVERY - TCRGO')}|
			+{'='*78}+
			{indent(__doc__, 3*'	')}
			"""
		),
		epilog=dedent(
			"""
			Developed by James Gatter & Sarah Nyquist of Shalek Lab (2020-2021). 
			Original Matlab pipeline written by Andy Tu and Tod Gierhann of Love Lab (2019).
			"""
		)
	)
	# REQUIRED ARGUMENTS
	parser.add_argument(
		"bam",
		metavar="SORTED_BAM",
		type=str,
		help="Path to sorted, single-end BAM"
	)
	parser.add_argument(
		'-f', "--fasta",
		type=str,
		required=True,
		help="Path to reference FASTA."
	)
	parser.add_argument(
		'-c', "--cdr3-positions-file",
		type=str,
		required=False,
		help=
			"Path to a two-column TSV containing reference names "
			"and the start V or end J base position. NOTE: If these "
			"positions are zero-indexed as opposed to one-indexed, "
			"also specify the --zero-indexed flag! If not provided, "
			"program attempts to identify the best ORF and positions."
	)
	parser.add_argument(
		'-z', "--zero-indexed",
		action="store_true",
		required=False,
		help=
			"If this and the CDR3 positions file are specified, "
			"the positions entered are treated as indices. "
			"and will not be reduced by 1."
	)
	# OPTIONAL ARGUMENTS
	parser.add_argument(
		'-q', "--query-list",
		type=str,
		required=False,
		help="The path to the query list file to read."
	)
	parser.add_argument(
		'-w', "--worker",
		type=int,
		metavar="#WORKER",
		default=1,
		help="The number ID of the worker which will read queries<W>.txt  (default: %(default)d)."
	)
	parser.add_argument(
		'-mf', "--minimum-frequency",
		type=float,
		default=0.0,
		help=
			"UMIs whose top VJ combination frequency is beneath the minimum frequency "
			"will not undergo CDR3 sequence reconstruction (default: %(default)d)."
	)
	parser.add_argument(
		'-mc', "--minimum-cdr3s",
		type=int,
		default=1,
		help=
			"UMIs with less than the minimum number of reads which have recovered CDR3s"
			"will not have their CDR3 statistics reported (default: %(default)d)."
	)
	parser.add_argument(
		'-r', "--exclude-relatives",
		action="store_true",
		help=
			"When identifying the VJ combination for a UMI, "
			"only take for CDR3 recovery reads whose top VJ is exactly equal "
			"to the highest counted segment variants. If not enabled, all reads "
			"are taken that match the highest counted segment's root (ex: TRAV24-1 "
			"is the top V for the UMI, but TRAV24-2 reads are also gathered for "
			"CDR3 recovery). Regardless the top segment's specific variant will be "
			"reported in the output information."
	)
	# TODO: Figure out why verbosity isn't working across modules.
	parser.add_argument(
		'-v', "--verbosity", 
		type=str,
		default="INFO",
		choices=["INFO", "VERBOSE", "WARNING", "SUCCESS", "ERROR"],
		help="The verbosity of the log. (default: %(default)s)."
	)
	parser.add_argument(
		'-i', "--input-path", 
		type=str,
		default="./out/queries/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	parser.add_argument(
		'-o', "--output-path", 
		type=str,
		default="./out/cdr3/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	args = parser.parse_args()
	main(args)