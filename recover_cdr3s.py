"""
This is the second script of the python pipeline. 
It is meant to take the filtered reads, find the top V and J alignments,
reconstruct the CDR3 sequences.
Run using python -m recover_cdr3s <args>

Requirements:
	Python 3.8.5, samtools, pysam 
"""
import argparse
import os.path as osp
from textwrap import dedent, indent
import pysam

from typing import List, Dict, Iterator, Set
from pathlib import Path

from tcrgo.bam import BAMDict, ReferenceDict
import tcrgo.io as io
from tcrgo import Log
log = Log(name=__name__)

def main(args):
	
	log.init(args.verbosity)
	log.info("Parsing input data...")
	bam = io.sort_and_index(args.bam, args.output_path)
	
	log.info("Reading CDR3 positions file and FASTA file...")
	cdr3_positions = None
	if args.cdr3_positions_file is not None:
		cdr3_positions = io.read_cdr3_positions(args.cdr3_positions_file)
	entries = io.parse_fasta(args.fasta)
	refdict = ReferenceDict()
	refdict.build(entries, cdr3_positions, args.zero_indexed)
	refdict.verify()

	log.info("Fetching filtered queries containing alignments to both V and J segments.")
	if args.query_list is not None:
		id_queries = io.read_query_list(args.query_list)
		args.worker = io.get_worker_id(args.query_list)
	else:
		id_queries = io.read_id_queries(args.input_path, args.worker)

	bamdict = BAMDict(bam)
	log.info("Finding the top V and J alignments for each query")
	# TODO: Consider splitting fully-mapped V/J reads and partially-mapped V/J reads into different camps
	# Perform de Bruijn graph sequence assembly on partially mapped reads.
	bamdict.build(id_queries, refdict)
	log.info("Finished retrieving top V and J alignments for each query.")
	log.info(f"Got {len(bamdict.get_barcodes())} barcode(s), "
			f"{len(bamdict.get_umis())} UMI(s) and "
			f"{len(bamdict.get_reads())} reads total")

	log.info(f"{len(bamdict.get_umis())}")
	for barcode in bamdict.get_barcodes():
		for umi in barcode.get_umis():
			umi.count_regions()
			umi.resolve_alignment_identities(refdict)
			reads_VJ = umi.get_top_VJ_reads(args.disclude_relatives)
			if len(reads_VJ) >= args.minimum_cdr3s \
				and umi.frequency_top_VJ >= args.minimum_frequency:
				for read in reads_VJ:
					read.get_cdr3()
				umi.select_cdr3(reads_VJ)

	log.info("Finished recovering CDR3 sequences for all VJ reads"
		 "and recovered top CDR3 for each UMI.")
	bamdict.write_cdr3_info(args.worker, args.output_path)
	log.info("Writing reference info for alignment tiebreaks.")
	bamdict.write_tiebreaks_alignments(args.worker, args.output_path)
	bamdict.close()
	log.success("DONE")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		prog="python -m recover_cdr3s", # Would be good to make a command-line alias and set that here
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
		type=str,
		help="Path to sorted, single-end BAM"
	)
	parser.add_argument(
		'-f', "--fasta",
		type=str,
		required=True,
		help="Path to reference FASTA."
	)
	parser.add_argument(
		'-c', "--cdr3-positions-file",
		type=str,
		required=False,
		help=
			"Path to a two-column TSV containing reference names "
			"and the start V or end J base position. NOTE: If these "
			"positions are zero-indexed as opposed to one-indexed, "
			"also specify the --zero-indexed flag!"
	)
	parser.add_argument(
		'-z', "--zero-indexed",
		action="store_true",
		required=False,
		help=
			"If this and the CDR3 positions file are specified, "
			"the positions entered are treated as indices. "
			"and will not be reduced by 1."
	)
	# OPTIONAL ARGUMENTS
	parser.add_argument(
		'-q', "--query-list",
		type=str,
		required=False,
		help="The path to the query list file to read."
	)
	parser.add_argument(
		'-w', "--worker",
		type=int,
		metavar="#WORKER",
		default=1,
		help="The number ID of the worker which will read queries<W>.txt  (default: %(default)d)."
	)
	parser.add_argument(
		'-mf', "--minimum-frequency",
		type=float,
		default=0.0,
		help=
			"UMIs whose top VJ combination frequency is beneath the minimum frequency "
			"will not undergo CDR3 sequence reconstruction (default: %(default)d)."
	)
	parser.add_argument(
		'-mc', "--minimum-cdr3s",
		type=int,
		default=1,
		help=
			"UMIs with less than the minimum number of reads which have recovered CDR3s"
			"will not have their CDR3 statistics reported (default: %(default)d)."
	)
	# TODO: This never happened. Perhaps repurpose for an additional check on the barcode level?
	# Ask Sarah if that's needed.
	parser.add_argument(
		'-s', "--second-pass",
		action="store_true",
		help=
			"After a top VJ combination is found, revisit non-top-VJ-combination reads "
			"to reclaim for CDR3 reconstruction reads which have primary OR high-scoring " 
			"secondary alignments to top V and J."
	)
	
	parser.add_argument(
		'-r', "--disclude-relatives",
		action="store_true",
		help=
			"When identifying the VJ combination for a UMI, "
			"only take for CDR3 recovery reads whose top VJ is exactly equal "
			"to the highest counted segment variants. If not enabled, all reads "
			"are taken that match the highest counted segment's root (ex: TRAV24-1 "
			"is the top V for the UMI, but TRAV24-2 reads are also gathered for "
			"CDR3 recovery). Regardless the top segment's specific variant will be "
			"reported in the output information."
	)
	# TODO: Figure out why verbosity isn't working across modules.
	parser.add_argument(
		'-v', "--verbosity", 
		type=str,
		default="INFO",
		choices=["INFO", "VERBOSE", "WARNING", "SUCCESS", "ERROR"],
		help="The verbosity of the log. (default: %(default)s)."
	)
	parser.add_argument(
		'-i', "--input-path", 
		type=str,
		default="./out/queries/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	parser.add_argument(
		'-o', "--output-path", 
		type=str,
		default="./out/cdr3/",
		help="The path to which the files from this program will output. (default: %(default)s )."
	)
	args = parser.parse_args()
	main(args)