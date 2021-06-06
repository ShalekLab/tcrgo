"""
Alignment and preprocessing via Drop-Seq Tools 2.4.0 and Bowtie2 version 2.4.1
Requires Samtools and pysam.

Note: `<basename>` is inputted using the `--basename BASENAME` argument.
	0. Input Data -> Unmapped BAM, starting from any of these substeps depending on the input:
		a. Raw BAM (`<basename>`.bam)
		b. Single-end FASTQ (`<basename>`.fastq)
		c. Pair-end FASTQs (`<basename>`_R1.fastq,`<basename>`_TCR.fastq)
	1. Unmapped BAM (`<basename>`_unmapped.bam) -> aligned and tagged BAM
		a. Tag cell barcodes (`<basename>_celltagged.bam`)
		b. Tag molecular barcodes (`<basename>_idtagged.bam`)
		c. Filter out reads with low quality cell/molecular barcodes (`<basename>`_filtered.bam)
		d. Trim 5’ primer sequence (`<basename>`_adaptertrimmed.bam)
		e. Trim 3’ polyA sequence (`<basename>`_trimmed.bam)
		f. SAM -> Fastq (`<basename>`_trimmed.fastq)
		g. Bowtie2 alignment (`<basename>`_aligned.sam)
		h. Sort STAR alignment in queryname order (`<basename>`_alignedsorted.bam)
		i. Merge STAR alignment tagged SAM to recover cell/molecular barcodes (`<basename>`_merged.bam)
		j. Add gene/exon and other annotation tags (`<basename>`_exontagged.bam)
		k. Barcode Repair
			i. Repair substitution errors (`<basename>_synthrepaired.bam`)
			ii. Repair indel errors (`<basename>`_repaired.bam)
			iii. Sort by coordinates (`<basename>`_repairedsorted.bam)
"""
import tcrgo.dropseq_tools as ds
import tcrgo.io as io
import argparse
from textwrap import dedent, indent
import os
import pysam

from tcrgo import Log
log = Log(name=__name__)

def main(args):
	log.init(args.verbosity)
	log.info("Verifying that Drop-Seq Tools, Picard, and Bowtie2 are all callable...")
	ds.test_tools(args.dropseq, args.picard, args.aligner)

	if not os.path.exists(args.output_path):
		os.makedirs(args.output_path)
	if args.basename is None:
		args.basename = os.path.basename(args.data[0].split('.')[0])

	# Files are produced in this order
	basename = os.path.join(args.output_path, args.basename)
	fastq_singleend = basename + ".fastq"
	fastq_barcode = basename + "_R1.fastq"
	fastq_biological = basename + "_TCR.fastq"
	bam_unmapped = basename + "_unmapped.bam"
	bam_idtagged = basename + "_idtagged.bam"
	bam_filtered = basename + "_filtered.bam"
	bam_trimmed = basename + "_trimmed.bam"
	fastq_trimmed = basename + "_trimmed.fastq"
	sam_aligned = basename + "_aligned.sam"
	bam_alignedsorted = basename + "_alignedsorted.bam"
	bam_merged = basename +"_merged.bam"
	bam_exontagged = basename + "_exontagged.bam"
	bam_synthrepaired = basename + "_synthrepaired.bam"
	bam_repaired = basename + "_repaired.bam"
	bam_repairedsorted = basename + "_repairedsorted.bam"

	log.info("Checking sequence data and FASTA arguments.") 
	if not os.path.isfile(args.fasta):
		log.error("Please enter a valid path to your FASTA file!")
	is_bam = io.is_bam(args.data)
	is_fastq_singleend = len(args.data) < 2
	if is_bam:
		bam_singleend = args.data[0]
	elif is_fastq_singleend:
		fastq_singleend = args.data[0]
	else:
		fastq_barcode, fastq_biological = args.data

	if is_bam and not os.path.isfile(fastq_singleend):
		log.info("Converting raw BAM to single-end FASTQ...")
		ds.bam2fq(bam_singleend, fastq_singleend)
	if is_fastq_singleend and not os.path.isfile(fastq_barcode) and not os.path.isfile(fastq_biological):
		log.info("Transforming Read1-index-1 FASTQ to Read1 and Read2 FASTQs")
		ds.transform_read_data(fastq_singleend, args.basename, fastq_barcode, fastq_biological, args.output_path)
	if not os.path.isfile(bam_unmapped):
		log.info("Converting R1 and R2 FASTQs to a single-end BAM")
		ds.fastq_to_bam(args.picard, fastq_barcode, fastq_biological, bam_unmapped, args.basename)
	if not os.path.isfile(bam_idtagged):
		log.info("Tagging BAM with cell barcode and UMI sequences...")
		ds.tag_identifiers_bam(args.dropseq, bam_unmapped, bam_idtagged)
	if args.filter_barcodes:
		if not os.path.isfile(bam_filtered):
			log.info("Filtering reads with low quality bases in cell barcode and UMI.")
			ds.filter_bam(args.dropseq, bam_idtagged, bam_filtered)
		if not os.path.isfile(bam_trimmed):
			log.info("Trimming adapter and poly A sequences...")
			#ds.trim_bam(args.dropseq, bam_idtagged, bam_trimmed)
			ds.trim_bam(args.dropseq, bam_filtered, bam_trimmed)
	else:
		if not os.path.isfile(bam_trimmed):
			log.info("Trimming adapter and poly A sequences...")
			ds.trim_bam(args.dropseq, bam_idtagged, bam_trimmed)
			#ds.trim_bam(args.dropseq, bam_filtered, bam_trimmed)

	if not os.path.isfile(sam_aligned):
		log.info(f"Aligning BAM using {args.aligner}")
		if args.aligner == "bowtie2":
			ds.bowtie2(bam_trimmed, args.fasta, sam_aligned)
		elif args.aligner == "star":
			if not os.path.isfile(fastq_trimmed):
				ds.bam2fq(bam_trimmed, fastq_trimmed)
			ds.star(fastq_trimmed, args.fasta, basename, sam_aligned)
	if not os.path.isfile(bam_alignedsorted):
		log.info("Sorting aligned SAM by queryname, outputting as BAM.")
		ds.sort_sam(args.picard, sam_aligned, bam_alignedsorted)
	if not os.path.isfile(bam_merged):
		log.info("Creating reference dictionary from FASTA for BAM merging.")
		ds.create_sequence_dictionary(args.picard, args.fasta)
		log.info("Merging the aligned BAM with the unmapped BAM to recover info for repair.")
		ds.merge_bam(args.picard, args.fasta, bam_trimmed, bam_alignedsorted, bam_merged)
	if not os.path.isfile(bam_exontagged):
		#log.info("Faux-tagging the merged BAM with dummy exon info.")
		log.info("Tagging the merged BAM with TRA/TRB exon info.")
		ds.fauxtag_exons(bam_merged, bam_exontagged)
	if not os.path.isfile(bam_repaired):
		log.info("Repairing insertion-deletion errors in the cell barcode sequences.")
		#ds.repair_bam(args.dropseq, bam_exontagged, bam_repaired, min_umis_per_cell=1)
		ds.bead_synthesis_errors(args.dropseq, bam_exontagged, bam_synthrepaired, min_umis_per_cell=1)
		log.info("Repairing substitution errors in the cell barcode sequences and collapsing barcodes.")
		ds.bead_substitution_errors(args.dropseq, bam_synthrepaired, bam_repaired, min_umis_per_cell=1)
		log.info("Sorting bam by coordinate, outputting final BAM.")
		pysam.sort("-o", bam_repairedsorted, bam_repaired)
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
		metavar="<BAM|FASTQ or FASTQ_R1 FASTQ_Index-1>",
		help="Path to single-end BAM/FASTQ or paths (delim. by comma) Read1-index-1 FASTQs"
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
		default="./out/alignment/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	parser.add_argument(
		'-a', "--aligner",
		type=str,
		default="bowtie2",
		choices=["bowtie2", "star"],
		required=False,
		help=
			"The aligner the program will choose for aligning the sequence data "
			"against the reference FASTA. (default: %(default)s )."
	)
	parser.add_argument(
		'-f', "--filter-barcodes",
		action="store_true",
		required=False,
		help=
			"Determines if the filter barcode based on quality step should be run. "
	)

	args = parser.parse_args()
	main(args)
