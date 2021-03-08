"""
This is the third and final script of the python pipeline. 
It is meant to aggregate the results from each recover_cdr3s worker
and summarize the information
Run using python -m summary <args>

Requirements:
	Python >3.8.5, samtools, pysam, biopython, pandas
"""
import argparse
from textwrap import dedent, indent
import os
import tcrgo.io as io
from tcrgo.collapse import collapse

from tcrgo import Log
log = Log('root')

def main(args):
	log.init(args.verbosity)
	if args.workers == "ALL":
		worker_range = io.list_cdr3_files(args.input_path)
	else:
		worker_range = [int(i) for i in args.workers.strip(':').split(':')]
		if worker_range[0] < 1:
			log.error("Please enter positive value for the start of the range.")
		if len(worker_range) > 1:
			worker_range = range(worker_range[0], worker_range[-1]+1)
		else:
			worker_range = range(1, worker_range[-1]+1)
	
	log.info("Aggregating CDR3 info files outputted from the recover_cdr3s script.")
	aggregated_cdr3_info = io.read_cdr3_info(worker_range, args.input_path)
	aggregated_cdr3_info = io.sort_and_index_dataframe(aggregated_cdr3_info, args.string_index)
	io.write_dataframe(aggregated_cdr3_info, args.output_path, "aggregated_cdr3_info.tsv")

	if args.collapse:
		log.info("Collapsing aggregated CDR3 info, this may take a long time.")
		outdir = os.path.join(args.output_path, "dgraphs")
		if args.plot:
			log.info(f"Plots for each topVJ combination will be plotted to {outdir}")
			if not os.path.isdir(outdir):
				os.makedirs(outdir)
		collapsed_cdr3_info = collapse(aggregated_cdr3_info, args.threshold, args.plot, outdir)
		collapsed_cdr3_info = io.sort_and_index_dataframe(collapsed_cdr3_info, args.string_index)
		io.write_dataframe(collapsed_cdr3_info, args.output_path, "aggrcollapsed_cdr3_info.tsv")

	log.info("Aggregating tie break information outputted from the recover_cdr3s script.")
	aggr_tiebreaks_alignments = io.read_tiebreaks_alignments(worker_range, args.input_path)
	io.write_dataframe(aggr_tiebreaks_alignments, args.output_path, "aggregated_tiebreaks_alignments.tsv")
	log.success("DONE")

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
		metavar="<'ALL' or NUMWORKERS or NUMFIRST:NUMLAST>",
		help=
			"The number or range [NUMFIRST, NUMLAST] of recover_cdr3s text "
			"files to be concatenated for summary. NUMFIRST must be > 0. ALL "
			"will automatically detect all files in the input directory."
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
		'-c', "--collapse", 
		action="store_true",
		help=
			"Perform BCUMI collapsing. BCUMIs are first binned by topVJ combination, then "
			"subgrouped by greater/equal to threshold (collapsers) vs. less than threshold "
			"(collapsees). For each collapser BCUMI, hamming distance comparisons are done "
			"with each collapsee BCUMI. If HD == 1, then CDR3 sequences of equal lengths "
			"are compared. If HD <= 1, then collapsees are collapsed into the collapser "
			"with the highest number of reads (or number of edges if there is a tie)."
			"aggrcollapsed_cdr3_info.tsv is created (default: %(default)s)."
	)
	parser.add_argument(
		'-p', "--plot", 
		action="store_true",
		help=
			"If BCUMI collapsing, plot Collapsers-Collapsees graphs for each reference combination."
			"If enabled, expect hundreds of graphs in 'output-path/dgraphs/ with cumulative size "
			"exceeding 100MB. The program will take longer to finish (default: %(default)s)."
	)
	parser.add_argument(
		'-t', "--threshold",
		type=int,
		default=10,
		help=
			"When BCUMI collapsing, BCUMIs with greater than or equal to this threshold will be " 
			"considered collapsers and the rest collapsees. Collapsees are merged into collapsers "
			"if the collapsing criteria is met. A lower number may result in more collapsing while "
  			"a higher number may be more strict and result in less collapsing (default: %(default)d)."
	)
	parser.add_argument(
		'-i', "--input-path",
		type=str,
		default="./out/cdr3/",
		help="Path to folder containing files outputted by recover_cdr3s.py (default: %(default)s)."
	)
	# TODO: change this
	parser.add_argument(
		'-o', "--output-path", 
		type=str,
		default="./out/cdr3/",
		help="The path to which the files from this program will output (default: %(default)s )."
	)
	args = parser.parse_args()
	main(args)