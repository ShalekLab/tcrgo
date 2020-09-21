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
import pysam

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
	
	log.info("Reading filtered queries containing V and J subregions")
	id_queries = io.read_id_queries(args.output_path, args.worker)

	bamdict = BAMDict(bam)
	log.info("Finding the top V and J alignments for each query")
	bamdict.build(id_queries)
	log.info("Finished retrieving top V and J alignments for each query.")

	log.info(f"Got {len(bamdict.get_barcodes())} barcode(s), "
			f"{len(bamdict.get_umis())} UMI(s) and "
			f"{len(bamdict.get_reads())} reads total")

	# TODO: Additional filtering?
	log.info(f"{len(bamdict.get_umis())}")
	for umi in bamdict.get_umis():
		umi.count_top_VJ()
		umi.count_top_VJ_by_alignments()

	log.info("Reading CDR3bases and FASTA files")
	cdr3_positions = io.read_cdr3_file(args.cdr3_positions_file)

	pysam.samtools.faidx(str(args.fasta))
	fasta = pysam.FastaFile(args.fasta)

	bamdict.reconstruct_cdr3s(fasta, cdr3_positions)
	
	log.info("Finished reconstructing CDR3 sequences for all reads")
	bamdict.cdr3_statistics()

	bamdict.close()
	log.success("Done")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m reconstruct_tcrs", # Would be good to make a command-line alias and set that here
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