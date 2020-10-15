"""
This is the third script of the python pipeline. 
It is meant to aggregate the results from each reconstruct_tcrs worker
and summarize the information
Run using python -m summary <args>

Requirements:
	Python 3.8.5, samtools, pysam, pandas
"""
import argparse
from pathlib import Path
import tcrgo.io as io
from textwrap import dedent, indent

from tcrgo import Log
log = Log(name=__name__)

def main(args):
	log.init(args.verbosity)
	worker_range = [int(i) for i in args.workers.strip(':').split(':')]
	if worker_range[0] < 1:
		log.error("Please enter positive value for the start of the range.")
	if len(worker_range) > 1:
		worker_range = range(worker_range[0], worker_range[-1]+1)
	else:
		worker_range = range(1, worker_range[-1]+1)

	aggregated_cdr3_info = io.read_cdr3_info(worker_range, args.input_path, args.string_index)
	io.write_aggregated_cdr3_info(aggregated_cdr3_info, args.output_path)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m summary", # Would be good to make a command-line alias and set that here
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
			Developed by James Gatter & Sarah Nyquist of Shalek Lab (2020). 
			Original Matlab pipeline written by Andy Tu and Tod Gierhann of Love Lab (2019).
			"""
		)
	)
	# REQUIRED ARGUMENTS
	# TODO: Perhaps make optional and autocollect all from input_directory
	parser.add_argument(
		"workers",
		type=str,
		metavar="<NUMWORKERS or NUMFIRST:NUMLAST>",
		help=
			"The number or range [NUMFIRST, NUMLAST] of reconstruct_tcrs text "
			"files to be concatenated for summary. NUMFIRST must be > 0."
	)
	# OPTIONAL ARGUMENTS
	parser.add_argument(
		'-v', "--verbosity", 
		type=str,
		default="INFO",
		choices=["INFO", "VERBOSE", "WARNING", "SUCCESS", "ERROR"],
		help="The verbosity of the log (default: %(default)s)."
	)
	parser.add_argument(
		'-s', "--string-index", 
		action="store_true",
		help=
			"Output the aggregated CDR3 info .tsv file with the BC and UMI strings"
			"as the indices instead of numbering them for easier reading (default: %(default)s)."
	)
	parser.add_argument(
		'-i', "--input-path",
		type=Path,
		default="./out/cdr3/",
		help="Path to folder containing files outputted by reconstruct_tcrs.py (default: %(default)s)."
	)
	# TODO: change this
	parser.add_argument(
		'-o', "--output-path", 
		type=Path,
		default="./out/cdr3/",
		help="The path to which the files from this program will output (default: %(default)s )."
	)
	args = parser.parse_args()
	main(args)