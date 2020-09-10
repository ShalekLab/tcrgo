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

def main(args):

	log.init(args.verbosity)
	log.info("Parsing input data...")
	bam = io.parse_data(args.data, args.output_path)
	log.info("Opened bam.")

	log.info("Fetching queries which contain alignments to V and J...")
	bamdict = BAMDict(bam)
	queries_VJ, queries_regions = bamdict.parse_queries()
	log.info(f"Fetched {len(queries_VJ)} reads.")
	bamdict.close()

	# TODO: What are the metrics by which we can filter down the list of queries?
	def filter_queries(queries_VJ, queries_regions):
		pass
	#filtered_queries = filter_queries(queries_VJ, queries_regions)
	#log.info("Finished filtering reads!")
	#io.output_queries(filtered_queries, args.output_path, args.divy_workers)
	log.info("Writing filtered queries containing V and J to file.")
	io.output_queries(queries_VJ, args.output_path, args.workers)
	
	log.success("Done!")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m filter_queries", # Would be good to make a command-line alias and set that here
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
		'-o', "--output-path", 
		type=Path,
		default="./tcrgo_out",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	parser.add_argument(
		'-w', "--workers",
		type=int,
		metavar="#Workers",
		default=1,
		help=
			"Divide the list of filtered reads by this number for processing "
			"by multiple HPCC tasks or VMs  (default: %(default)d)."
	)
	args = parser.parse_args()
	main(args)