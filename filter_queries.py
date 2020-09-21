"""
This is the first script of the python pipeline. 
It is meant to filter reads down to only those which have V and J alignments.
Run using python -m filter_queries <args>

Requirements:
	Python 3.8.5, samtools, pysam 
"""
import argparse
import os.path as osp
from textwrap import dedent, indent

from typing import List, Dict, Iterator, Set
from pathlib import Path

from tcrgo.bam import BAMDict
import tcrgo.io as io
from tcrgo import Log
log = Log(name=__name__)

# TODO: Seems like there ARE standard tags for cell barcode/UMI. CR and OX/MI

def main(args):

	log.init(args.verbosity)
	log.info("Parsing input data...")
	bam = io.parse_data(args.data, args.output_path)
	log.info("Fetching unique queries which contain alignments to V and J and grouping them by Barcode-UMI...")
	queries_VJ, queries_nonVJ = io.parse_queries(bam)
	log.info(f"Fetched {len(queries_VJ)} VJ reads and {len(queries_nonVJ)} non-VJ reads.")
	bam.close()

	# TODO: Require user do reference sequence similarity checks for subregion variants.
	# 		We will assume for now dissimilar subregion variants.
	#		It would be useful to report pairs of regions that frequently appear together
	#		as alignments of a read.
	#		Maybe could use pysam.PileUpColumn/PileUp for this?
	
	log.info("Writing all queries not containing V and J alignments to file.")
	io.output_nonVJ(queries_nonVJ, args.output_path)
	log.info(f"Filtering out UMIs with VJ below {args.minimum_reads} reads and randomly "
			f"subsampling those with reads exceeding {args.maximum_reads} reads.")
	log.info(f"Will be writing filtered queries containing V and J to {args.workers} file(s).")
	id_queries = io.filter_queries(queries_VJ, args.minimum_reads, args.maximum_reads) # Split by delim '|' to convert to dictionary of BC|UMI : QNAME
	io.output_grouped_VJ_queries(id_queries, args.workers, args.output_path)
	log.success("Done!")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m filter_queries", # Would be good to make a command-line alias and set that here
		usage="",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=dedent(
			f"""
			+{'='*78}+
			|{'{:^78}'.format('SEQ-WELL TCR RECOVERY - TCRGO?')}|
			+{'='*78}+
			{indent(__doc__, 3*'	')}
			"""
		),
		epilog=dedent(
			"""
			Developed by James Gatter & Sarah Nyquist of Shalek Lab (2020). 
			Original Matlab pipeline written by Andy Tu and Tod Gierhann of Love Lab (2019).
			"""
		)
	)
	# REQUIRED ARGUMENTS
	parser.add_argument(
		"data", 
		type=Path,
		nargs='+',
		metavar="<BAM or FASTQ_R1 FASTQ_R2>",
		help="Path to single-end BAM or the paired-end FASTQs"
	)

	# OPTIONAL ARGUMENTS
	parser.add_argument(
		'-v', "--verbosity", 
		type=str,
		default="INFO",
		choices=["INFO", "VERBOSE", "WARNING", "SUCCESS", "ERROR"],
		help="The verbosity of the log. (default: %(default)s)."
	)
	parser.add_argument(
		'-lt', "--minimum_reads",
		type=int,
		default=5,
		help="UMIs with less than the minimum number of reads will be filtered out (default: %(default)d)."
	)
	parser.add_argument(
		'-gt', "--maximum_reads",
		type=int,
		default=1000,
		help=
			"UMIs with more than the maximum number of reads will have their reads"
			"randomly subsampled down to the maximum (default: %(default)d)."
	)
	parser.add_argument(
		'-o', "--output-path", 
		type=Path,
		default="./tcrgo_out",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	parser.add_argument(
		'-w', "--workers",
		type=int,
		metavar="#WORKERS",
		default=1,
		help=
			"Divide the list of filtered reads by this number for processing "
			"by multiple HPCC tasks or VMs  (default: %(default)d)."
	)
	args = parser.parse_args()
	main(args)