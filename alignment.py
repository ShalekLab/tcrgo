"""Alignment and preprocessing via Drop-Seq Tools 2.4.0 and Bowtie2 version 2.4.1"""
import tcrgo.dropseq_tools as ds
import argparse
from textwrap import dedent, indent
from pathlib import Path
import os

from tcrgo import Log
log = Log(name=__name__)

def main(args):
	"""
	1. Unmapped BAM -> aligned and tagged BAM
		a. Tag cell barcodes
		b. Tag molecular barcodes
		c. Trim 5’ primer sequence
		d. Trim 3’ polyA sequence
		e. SAM -> Fastq
		f. STAR alignment
		g. Sort STAR alignment in queryname order
		h. Merge STAR alignment tagged SAM to recover cell/molecular barcodes
		i. Add gene/exon and other annotation tags
		j. Barcode Repair
			i. Repair substitution errors (DetectBeadSubstitutionErrors)
			ii. Repair indel errors (DetectBeadSynthesisErrors)
	"""

	log.info("Verifying that Drop-Seq Tools, Picard, and Bowtie2 are all callable...")
	ds.test_tools(args.dropseq, args.picard)

	log.info("Checking sequence data and FASTA arguments.") 
	if not os.path.isfile(args.fasta):
		log.error("Please enter a valid path to your FASTA file!")
	"""
	if isinstance(args.data, list):
		if len(args.data) != 2:
			log.error("For input sequence data, please enter either an unmapped BAM "
				f"or paired end FASTQs; FASTQ_R1, FASTQ_R2. Received {len(args.data)}.")
		log.info(f"Converting {args.data[0]} and {args.data[1]} to BAM")
		bam_unmapped = ds.fastq_to_bam(args.picard, args.data[0], args.data[1], bam_unmapped)
	else:
		bam_unmapped = args.data
	"""
	if len(args.data) > 2:
		log.error("For input sequence data, please enter either an unmapped BAM "
			f"or paired end FASTQs; FASTQ_R1, FASTQ_R2. Received {len(args.data)}.")
	elif len(args.data) == 2:
		log.info(f"Converting {args.data[0]} and {args.data[1]} to BAM")
		bam_unmapped = ds.fastq_to_bam(args.picard, args.data[0], args.data[1], bam_unmapped)
	else:
		bam_unmapped = args.data[0]

	if not os.path.exists(args.output_path):
		os.makedirs(args.output_path)
	if args.basename is None:
		args.basename = os.path.basename(bam_unmapped.split('.')[0])
	basename = os.path.join(args.output_path, args.basename)
	bam_idtagged = basename + "_idtagged.bam"
	bam_trimmed = basename + "_trimmed.bam"
	sam_aligned = basename + "_aligned.sam"
	bam_alignedsorted = basename + "_alignedsorted.bam"
	bam_merged = basename +"_merged.bam"
	bam_exontagged = basename + "_exontagged.bam"
	bam_repaired = basename + "_repaired.bam"

	log.info("Tagging BAM with cell barcode and UMI sequences...")
	ds.tag_identifiers_bam(args.dropseq, bam_unmapped, bam_idtagged)
	log.info("Trimming adapter and poly A sequences...")
	ds.trim_bam(args.dropseq, bam_idtagged, bam_trimmed)
	log.info("Aligning BAM using Bowtie2")
	ds.align_bam(bam_trimmed, args.fasta, sam_aligned)
	log.info("Sorting aligned SAM by queryname, outputting as BAM.")
	ds.sort_sam(args.picard, sam_aligned, bam_alignedsorted)
	log.info("Creating reference dictionary from FASTA for BAM merging.")
	ds.create_sequence_dictionary(args.picard, args.fasta)
	log.info("Merging the aligned BAM with the unmapped BAM to recover info for repair.")
	ds.merge_bam(args.picard, args.fasta, bam_trimmed, bam_alignedsorted, bam_merged)
	log.info("Faux-tagging the merged BAM with dummy exon info.")
	ds.fauxtag_exons(bam_merged, bam_exontagged)
	log.info("Repairing substitution and indel errors in the BAM barcode sequences.")
	ds.repair_bam(args.dropseq, bam_exontagged, bam_repaired, min_umis_per_cell=1)
	log.success("DONE")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m alignment", # Would be good to make a command-line alias and set that here
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
		type=str,
		nargs='+',
		metavar="<BAM or FASTQ_R1 FASTQ_R2>",
		help="Path to single-end BAM or the paired-end FASTQs"
	)
	parser.add_argument(
		'-d', "--dropseq", 
		type=str,
		help="Path to the Drop-Seq Tools jar file. Recommended you use a version 2.3.0 or 2.4.0."
	)
	parser.add_argument(
		'-p', "--picard", 
		type=str,
		help=
			"Path to the Picard jar file. Recommended you use a version of Picard that is "
			"bundled with/compatible with your version of Drop-Seq Tools."
	)
	parser.add_argument(
		'-f', "--fasta",
		type=str,
		required=True,
		help="Path to reference FASTA."
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
		'-b', "--basename",
		type=str,
		default=None,
		help=
			"The basename to affix to all output files. If not specified, "
			"Python's os module will be used to ascertain."
	)
	parser.add_argument(
		'-o', "--output-path", 
		type=str,
		default="./tcrgo_out/alignment/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	args = parser.parse_args()
	main(args)