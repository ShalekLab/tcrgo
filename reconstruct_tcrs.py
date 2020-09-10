"""
This is the second script of the python pipeline. 
It is meant to take the filtered reads, find the top V and J alignments,
reconstruct the CDR3 sequences.
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
	bam = io.sort_and_index(args.bam, args.output_path)
	log.info("Opened bam.")
	
	log.info("Fetching filtered queries containing V and J subregions")
	queries = io.input_queries(args.output_path, args.worker)

	bamdict = BAMDict(bam)
	log.info("Finding the top V and J alignments for each query")
	bamdict.parse_bam(queries)
	log.info("Finished retrieving top V and J alignments for each query.")

	# TODO: Additional filtering?

	log.info("Reading CDR3bases and FASTA files")
	cdr3_positions = io.read_cdr3_file(args.cdr3_positions_file)

	pysam.samtools.faidx(args.fasta)
	fasta = pysam.FastaFile(args.fasta)

	for read in bamdict.get_reads():
		CDR3_end = read.get_cdr3_position(read.top_J, cdr3_positions)
		CDR3_start = read.get_cdr3_position(read.top_V, cdr3_positions)

		# TODO: J and V ref start and end actually needed here? Or do we want the whole?
		seq_ref_J = fasta.fetch(
			read.top_J.reference_name,
			read.top_J.reference_start,
			read.top_J.reference_end
		)
		seq_ref_V = fasta.fetch(
			read.top_V.reference_name,
			read.top_V.reference_start,
			read.top_V.reference_end
		)

		if read.top_J.query_alignment_end - CDR3_end > 0:
			sequence = seq_ref_V + read.top_V.
		else:

	bamdict.close()
	log.success("Done")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m process_barcodes", # Would be good to make a command-line alias and set that here
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
		"bam",
		metavar="SORTED_BAM",
		type=Path,
		help="Path to sorted, single-end BAM"
	)
	parser.add_argument(
		'-f', "--fasta",
		type=Path,
		required=True,
		help="Path to reference FASTA."
	)
	parser.add_argument(
		'-c', "--cdr3-positions-file",
		type=Path,
		required=True,
		help=
			"Path to a two-column TSV containing reference names "
			"and the start V or end J base position"
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
		'-w', "--worker",
		type=int,
		metavar="#WORKER",
		default=1,
		help="The number ID of the worker which will read queries<W>.txt  (default: %(default)d)."
	)
	args = parser.parse_args()
	main(args)