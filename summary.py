"""
This is the third script of the python pipeline. 
It is meant to aggregate the results from each reconstruct_tcrs worker
and summarize the information
Run using python -m summary <args>

Requirements:
	Python 3.8.5, samtools, pysam
"""
import argparse

def main(args):
file.write(
	"BC_index\tUMI_index\tBC\tUMI\tnReads\t"
	"topV_region\ttopV_nReads\ttopV_freq\t"
	"topJ_region\ttopJ_nReads\ttopJ_freq\t
	"CDR3_nuc\tCDR3_nReads\tCDR3_freq\t"
	"TRAV_nReads\tTRAJ_nReads\tTRAC_nReads\t"
	"TRBV_nReads\tTRBJ_nReads\tTRBC_nReads\t\n"
)

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
	parser.add_argument(
		"input_directory",
		metavar="INPUT_DIRECTORY",
		type=Path,
		help="Path to folder containing files outputted by reconstruct_tcrs."
	)
	# OPTIONAL ARGUMENTS
	parser.add_argument(
		'-v', "--verbosity", 
		type=str,
		default="INFO",
		choices=["INFO", "VERBOSE", "WARNING", "SUCCESS", "ERROR"],
		help="The verbosity of the log. (default: %(default)s)."
	)
	# TODO: change this
	parser.add_argument(
		'-o', "--output-path", 
		type=Path,
		default="./tcrgo_out",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	# TODO: Handle this later
	parser.add_argument(
		'-w', "--workers",
		type=str,
		metavar="NUMFIRST:NUMLAST",
		default="ALL",
		help=\
			"The range of reconstruct_tcrs text files to be summarized summary (default: %(default)d)."
			"If not entered, all text files in INPUT_DIRECTORY will be collected for summary."
	)
	args = parser.parse_args()
	main(args)