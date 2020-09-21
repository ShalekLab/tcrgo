"""A hack for statistical analysis"""

import argparse
import os.path as osp
from textwrap import dedent, indent
import pysam

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

	cdr3_positions = io.read_cdr3_file(args.cdr3_positions_file)

	pysam.samtools.faidx(str(args.fasta))
	fasta = pysam.FastaFile(args.fasta)

	log.info("Parsing all queries into BAM Dict")
	bamdict = BAMDict(bam)
	bamdict.parse_queries_full()

	log.info("Writing all unique queries to file.")
	bamdict.write_reads(fasta, cdr3_positions)
	
	fasta.close()
	bamdict.close()
	log.success("Done!")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m statistics", # Would be good to make a command-line alias and set that here
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