"""
Subprocess calls to Drop-seq Tools
"""
import subprocess as sp
import pysam
import os

from pathlib import Path
from typing import List, Tuple, Dict
from textwrap import dedent

from tcrgo.log import Log
log = Log(name=__name__)
log.formatter.datefmt = "%H:%M:%S"

def execute(command: str): #-> str:
	try:
		output = sp.check_output(command.split())
	except sp.CalledProcessError:
		log.error(
			f"{command.split()} terminated with non-zero exit status, "
			"please check the above messages for more detail."
		)
	return output

def test_tools(dropseq_jar: str, picard_jar: str):
	"""Notify user if tools could not be called from commandline"""
	#execute(f"java -jar {dropseq_jar} TagBamWithReadSequenceExtended --version")
	#execute(f"java -jar {picard_jar} SortSam --version")
	if not os.path.isfile(dropseq_jar):
		log.error(f"Could not find {dropseq_jar}")
	if not os.path.isfile(picard_jar):
		log.error(f"Could not find {picard_jar}")
	execute("bowtie2 --version")

def fastq_to_bam(picard_jar: str, fastq_barcode: str, fastq_biological: str, bam_unmapped: str) -> str:
	"""Convert FASTQ R1 and R2 to unmapped BAM using Picard's FastqToSam"""
	command = \
		f"""
		java -Dsamjdk.compression_level=6 -Xmx4g -jar {picard_jar} FastqToSam \
			OUTPUT={bam_unmapped} \
			USE_SEQUENTIAL_FASTQS=true \
			FASTQ={fastq_barcode} \
			FASTQ2={fastq_biological} \
			QUALITY_FORMAT=Standard \
			SORT_ORDER=queryname
		""" # Might need "SAMPLE_NAME" arg?
	log.sep()
	execute(command)
	log.sep()
	return bam_unmapped

def tag_identifiers_bam(dropseq_jar: str, bam_untagged: str, bam_idtagged: str) -> str:
	"""
	Tag the BAM with tags for the UMI + Barcode
	bam_untagged = The unprocessed, unmapped BAM.
	bam_tagged = the cell- and UMI-tagged BAM.
	"""
	dirname = os.path.dirname(bam_idtagged)
	bam_celltagged = os.path.basename(bam_idtagged).replace("_idtagged.bam", "_celltagged.bam")
	bam_celltagged = os.path.join(dirname, bam_celltagged)
	summary_tagged_cellular = os.path.join(dirname, "summary_tagged_cellular.txt")
	summary_tagged_molecular = os.path.join(dirname, "summary_tagged_molecular.txt")
	command = \
		f"""
		java -Dsamjdk.compression_level=1 -Xmx4000m -jar {dropseq_jar} \
		TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_untagged} \
			OUTPUT={bam_celltagged} \
			BARCODE_QUALITY_TAG=CY \
			SUMMARY={summary_tagged_cellular} \
			BASE_RANGE="1-12" \
			BASE_QUALITY=10 \
			BARCODED_READ=1 \
			DISCARD_READ=false \
			TAG_NAME=XC \
			NUM_BASES_BELOW_QUALITY=1
		"""
	log.sep()
	execute(command)
	log.sep()
	command = \
		f"""
		java -Dsamjdk.compression_level=6 -Xmx4000m -jar {dropseq_jar} \
			TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_celltagged} \
			OUTPUT={bam_idtagged} \
			BARCODE_QUALITY_TAG=UY \
			SUMMARY={summary_tagged_molecular} \
			BASE_RANGE="13-20" \
			BASE_QUALITY=10 \
			BARCODED_READ=1 \
			DISCARD_READ=true \
			TAG_NAME=XU \
			NUM_BASES_BELOW_QUALITY=1
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_idtagged

def trim_bam(dropseq_jar: str, bam_idtagged: str, bam_trimmed: str) -> str:
	if not bam_idtagged.endswith("_idtagged.bam"):
		log.error("Ensure bam_tagged ends in '_idtagged.bam'")
	#bam_adaptertrimmed = bam_idtagged.replace("_idtagged.bam", "_adaptertrimmed.bam")
	
	dirname = os.path.dirname(bam_trimmed)
	bam_adaptertrimmed = os.path.basename(bam_idtagged).replace("_idtagged.bam", "_adaptertrimmed.bam")
	bam_adaptertrimmed = os.path.join(dirname, bam_adaptertrimmed)
	summary_adaptertrimmed = os.path.join(dirname, "summary_trimming_adapter.txt")
	summary_polyAtrimmed = os.path.join(dirname, "summary_trimming_polyA.txt")
		
	command = \
		f"""
		java -Dsamjdk.compression_level=1 -Xmx4g -jar {dropseq_jar} 
		TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_idtagged} \
			OUTPUT={bam_adaptertrimmed} \
			OUTPUT_SUMMARY={summary_adaptertrimmed} \
			SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
			MISMATCHES=0 \
			NUM_BASES=5
		"""
	log.sep()
	execute(command)
	log.sep()
	command = \
		f"""
		java -Dsamjdk.compression_level=2 -Xmx4g -jar {dropseq_jar}
		PolyATrimmer VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_adaptertrimmed} \
			OUTPUT={bam_trimmed} \
			OUTPUT_SUMMARY={summary_polyAtrimmed} \
			ADAPTER=^XU^XCACGTACTCTGCGTTGCTACCACTG \
			MISMATCHES=0 \
			NUM_BASES=6 \
			USE_NEW_TRIMMER=true
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_trimmed

@log.time
def remove_scoretags(bam: str) -> str:
	"""
	Removes any empty AS and XS tags prior to alignment so only 
	cell and UMI tags remain after bowtie2 --preserve-tags
	"""
	log.info("Removing stray score tags from BAM to prevent duplication in alignment output.")
	pysam.index(bam)
	bam_detagged = bam.replace("_trimmed.bam", "_detagged.bam")
	if bam_detagged == bam:
		log.error("Ensure the input bam ends with '_trimmed.bam'")
	with pysam.AlignmentFile(bam, 'rb') as tagged:
		with pysam.AlignmentFile(bam_detagged, 'wb', template=tagged) as detagged:
			for query in tagged:
				tags = list()
				for tag, value in query.get_tags():
					if tag == "XS" or tag == "AS":
						continue
					tags.append((tag, value))
				query.tags = tuple(tags)
				detagged.write(query)
	return bam_detagged
	

def align_bam(bam_trimmed: str, fasta: str, sam_aligned: str) -> str:
	"""Convert SAM to FASTQ using Picard then align to reference using STAR"""
	if not sam_aligned.endswith("_aligned.sam"):
		log.error("sam_aligned must end in '_aligned.sam'.")
	bam_detagged = remove_scoretags(bam_trimmed) # Avoid preserving empty score tags
	bam_sorted = bam_detagged.replace("_detagged.bam", "_sorted.bam")
	pysam.sort("-n", "-o", bam_sorted, bam_detagged) # Sort by QNAME
	pysam.index(bam_sorted)

	index_prefix = os.path.basename(fasta.split('.')[0])
	index_prefix = os.path.join(os.path.dirname(sam_aligned), index_prefix)
	command = f"bowtie2-build -f {fasta} --threads 32 {index_prefix}"
	log.sep()
	execute(command)
	command = f"bowtie2 -x {index_prefix} -a --very-sensitive-local -p 32 -S {sam_aligned} -b {bam_sorted} --preserve-tags"
	execute(command)
	log.sep()
	return sam_aligned

def sort_sam(picard_jar: str, sam_aligned: str, bam_alignedsorted: str) -> str:
	bam_alignedsorted = sam_aligned.replace("_aligned.sam", "_alignedsorted.bam")
	command = \
		f"""
		java -Dsamjdk.compression_level=1 -Xms4000m -jar {picard_jar} \
		SortSam VALIDATION_STRINGENCY=SILENT \
			INPUT={sam_aligned} \
			OUTPUT={bam_alignedsorted} \
			MAX_RECORDS_IN_RAM=2000000 \
			SORT_ORDER=queryname
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_alignedsorted

def create_sequence_dictionary(picard_jar: str, fasta: str) -> str:
	if fasta.endswith(".fa"):
		fasta_dict = fasta.replace(".fa", '.dict')
	elif fasta.endswith(".fasta"):
		fasta_dict = fasta.replace(".fasta", '.dict')
	if os.path.isfile(fasta_dict):
		log.warn(f"{fasta_dict} already exists!")
	else:
		command = \
			f"""
			java -Xms4000m -jar {picard_jar} \
			CreateSequenceDictionary \
				R={fasta} \
				O={fasta_dict}
			"""
		log.sep()
		execute(command)
		execute(f"chmod a+rw {fasta_dict}")
		log.sep()
	return fasta_dict

def merge_bam(picard_jar: str, fasta: str, bam_trimmed: str, bam_alignedsorted: str, merged_bam: str):
	merged_bam = bam_trimmed.replace("_trimmed.bam", "_merged.bam")
	command = \
		f"""
		java -Dsamjdk.compression_level=1 -Xms4000m -jar {picard_jar} \
			MergeBamAlignment VALIDATION_STRINGENCY=SILENT \
			ALIGNED_BAM={bam_alignedsorted} \
			UNMAPPED_BAM={bam_trimmed} \
			OUTPUT={merged_bam} \
			REFERENCE_SEQUENCE={fasta} \
			ATTRIBUTES_TO_RETAIN=XS \
			INCLUDE_SECONDARY_ALIGNMENTS=true \
			PAIRED_RUN=false \
			SORT_ORDER=coordinate \
			CLIP_ADAPTERS=false
		""" # There was a q after the last 'false'?
	log.sep()
	execute(command)
	log.sep()
	return merged_bam

def fauxtag_exons(bam_merged: str, bam_exontagged: str) -> str:
	log.sep()
	pysam.index(bam_merged)
	log.info(f"Faux-tagging {bam_merged} to produce {bam_exontagged}.")
	with pysam.AlignmentFile(bam_merged, 'rb') as untagged:
		with pysam.AlignmentFile(bam_exontagged, 'wb', template=untagged) as tagged:
			exon_tags = [("gf", "CODING", "Z"), ("gn", "GENEA", "Z"), ("gs", "+", "Z")]
			for query in untagged:
				tags = list(query.get_tags())
				query.set_tags(tags+exon_tags)
				tagged.write(query)
	log.sep()
	return bam_exontagged

def repair_bam(dropseq_jar: str, bam_exontagged: str, bam_repaired: str, min_umis_per_cell: int=1):
	"""
	Identifies cell barcodes with
	aberrant "fixed" UMI bases. If only the last UMI base is fixed as a T, the cell barcode is corrected (the
	last base is trimmed off) and all cell barcodes with identical sequence at the first 11 bases are merged
	together. If any other UMI base is fixed, the reads with that cell barcode are discarded.
	"""
	if not bam_exontagged.endswith("_exontagged.bam"):
		log.error("Exon-tagged BAM must end with '_exontagged.bam'")
	dirname = os.path.dirname(bam_repaired)
	bam_subrepaired = os.path.basename(bam_exontagged).replace("_exontagged.bam", "_subrepaired.bam")
	bam_subrepaired = os.path.join(dirname, bam_subrepaired)
	report_substitution_error = os.path.join(dirname, "report_substitution_error.txt")
	stats_synthesis_error = os.path.join(dirname, "stats_synthesis_error.txt")
	summary_synthesis_error = os.path.join(dirname, "summary_synthesis_error.txt")
	report_synthesis_error = os.path.join(dirname, "report_synthesis_error.txt")
	command = \
		f"""
		java -Dsamjdk.compression_level=2 -Xmx4g -jar {dropseq_jar} \
		DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
			VERBOSITY=DEBUG \
			INPUT={bam_exontagged} \
			OUTPUT={bam_subrepaired} \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			MOLECULAR_BARCODE_TAG=XU \
			NUM_THREADS=4 \
			OUTPUT_REPORT={report_substitution_error}
		"""
	log.sep()
	execute(command)
	log.sep()
	command = \
		f"""
		java -Dsamjdk.compression_level=1 -Xmx4g -jar {dropseq_jar} \
		DetectBeadSynthesisErrors VALIDATION_STRINGENCY=SILENT \
			VERBOSITY=DEBUG \
			INPUT={bam_subrepaired} \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			OUTPUT_STATS={stats_synthesis_error} \
			SUMMARY={summary_synthesis_error} \
			REPORT={report_synthesis_error} \
			MOLECULAR_BARCODE_TAG=XU \
			NUM_THREADS=4 \
			OUTPUT={bam_repaired}
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_repaired

# TODO: Python untested.
def collapsebarcodesinplace(dropseq_jar: str, bam: str, bam_output: str):
	command = \
		f"""
		java -jar {dropseq_jar} CollapseBarcodesInPlace \
			INPUT={bam} \
			OUTPUT={bam_output} \
			PRIMARY_BARCODE=XC \
			EDIT_DISTANCE=1 \
			FIND_INDELS=true \
			OUT_BARCODE=XC \
			MINIMUM_MAPPING_QUALITY=10 \
			MIN_NUM_READS_CORE=100 \
			MIN_NUM_READS_NONCORE=1 \
			FILTER_PCR_DUPLICATES=false \
			USE_JDK_INFLATER=true \
			USE_JDK_DEFLATER=true \
			NUM_THREADS=4
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_output

# TODO: Python untested.
def collapsetagwithcontext(dropseq_jar: str, bam: str, bam_output: str):
	command = \
		f"""
		java -Xmx4g -jar {dropseq_jar} CollapseTagWithContext \
			VERBOSITY=DEBUG \
			INPUT={bam} \
			COLLAPSE_TAG=XU \
			CONTEXT_TAGS=XC \
			OUT_TAG=XU \
			OUTPUT={bam_output} \
			EDIT_DISTANCE=1 \
			FIND_INDELS=true \
			READ_MQ=10
			MIN_COUNT=1 \
			DROP_SMALL_COUNTS=false \
			NUM_THREADS=4 \
			USE_JDK_INFLATER=true \
			USE_JDK_DEFLATER=true
		""" # COUNT_TAGS option?
	log.sep()
	execute(command)
	log.sep()
	return bam_output

