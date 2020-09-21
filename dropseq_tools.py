"""
Subprocess calls to Drop-seq Tools
"""
import subprocess as sp

from pathlib import Path
from typing import List, Tuple, Dict
#import pysam
#BAM = pysam.libcalignmentfile.AlignmentFile

from log import Log
log = Log(name=__name__)
log.formatter.datefmt = "%H:%M:%S"

# TODO
def test_path():
	"""Notify user if tools could not be called from commandline"""
	pass


def fastq_to_bam(fastq_barcodes: Path, fastq_biological: Path) -> Path:
	"""Convert FASTQ R1 and R2 to unmapped BAM using Picard's FastqToSam"""
	
	#if str(fastq_barcodes).endswith(".fastq.gz"):
	#	bam = fastq_barcodes.replace(".fastq.gz", ".bam")
	#elif str(fastq_barcodes).endswith(".fastq"):
	#	bam = fastq_barcodes.replace(".fastq", ".bam")
	#else:
	#	log.error("FASTQ does not end with either '.fastq' or '.fastq.gz'")
	
	log.info("Running Picard's FastqtoSam to produce an unmapped BAM.")
	command = [
		"java", "-Dsamjdk.compression_level=6", "-Xmx4G", "-jar", "picard.jar", "FastqToSam",
			f"OUTPUT={bam}",
			"USE_SEQUENTIAL_FASTQS=true",
			f"FASTQ={fastq_barcodes}",
			f"FASTQ2={fastq_biological}",
			"QUALITY_FORMAT=Standard",
			#"SAMPLE_NAME=~{sample_id}",
			"SORT_ORDER=queryname"
			#">", "FastqToSam.log"
	]
	log.verbose(' '.join(command))
	try:
		sp.check_call(command)
	except sp.CalledProcessError:
		log.error("Picard FastqToSam terminated with non-zero exit status, "
				"please check the FastqToSam.log for more details.")
	return bam

def add_tags_barcodes(bam: Path) -> Path:
	"""Tag the BAM with tags for the UMI + Barcode"""
	"""
	java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=~{unmapped_bam} \
		OUTPUT=pipe1 \
		~{true="BARCODE_QUALITY_TAG=CY" false="" quality_tags} \
		SUMMARY=~{sample_id}_tagged_cellular_summary.txt \
		BASE_RANGE=~{cellular_barcode_base_range} \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=false \
		TAG_NAME=XC \
		NUM_BASES_BELOW_QUALITY=1 | \
		
	java -Dsamjdk.compression_level=6 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe1 \
		OUTPUT="~{sample_id}_tagged.bam" \
		~{true="BARCODE_QUALITY_TAG=UY" false="" quality_tags} \
		SUMMARY="~{sample_id}_tagged_molecular_summary.txt" \
		BASE_RANGE=~{umi_base_range} \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=true \
		TAG_NAME=XM \
		NUM_BASES_BELOW_QUALITY=1
	"""
	pass

def filter_bam(bam:Path) -> Path:
	"""
	"""
	"""
	java -Dsamjdk.compression_level=1 -Xmx3000m -jar /software/Drop-seq_tools/jar/dropseq.jar FilterBam VALIDATION_STRINGENCY=SILENT \
		INPUT=~{unmapped_bam} \
		OUTPUT="~{sample_id}_filtered.bam" \
		TAG_REJECT=XQ
	"""
	pass

def trim_bam(bam: Path) -> Path:
	"""
	"""
	"""
	java -Dsamjdk.compression_level=1 -Xmx$command_mem_mb -jar /software/Drop-seq_tools/jar/dropseq.jar TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
		INPUT=~{unmapped_bam} \
		OUTPUT=pipe1 \
		OUTPUT_SUMMARY="~{sample_id}_adapter_trimming_report.txt" \
		SEQUENCE=~{trim_sequence} \
		MISMATCHES=0 \
		NUM_BASES=~{trim_num_bases} | \
		java -Dsamjdk.compression_level=2 -Xmx$command_mem_mb -jar /software/Drop-seq_tools/jar/dropseq.jar PolyATrimmer VALIDATION_STRINGENCY=SILENT \
		INPUT=pipe1 \
		OUTPUT="~{sample_id}_trimmed.bam" \
		OUTPUT_SUMMARY="~{sample_id}_polyA_trimming_report.txt" \
		MISMATCHES=0 \
		NUM_BASES=6 \
		NEW=true
	"""
	pass

def star(fastq: Path) -> Path:
	"""Convert SAM to FASTQ using Picard then align to reference using STAR"""
	"""java -Xmx5000m -jar /software/picard.jar SamToFastq \
      VALIDATION_STRINGENCY=SILENT \
      INPUT=~{input_bam} \
      FASTQ=/dev/stdout | \
	STAR \
      --genomeDir genome_dir \
      --outStd Log \
      --readFilesIn /dev/stdin \
      --runThreadN ~{cpu} \
      --outSAMtype BAM Unsorted \
      --outFileNamePrefix "~{sample_id}_" \
      ~{" " + flags}
	"""
	pass

def add_tags_genes_exons(aligned_bam: Path, unmapped_bam: Path) -> Path:
	"""Sort and Merge BAM Add gene/exon and other annotation tags"""
	"""
	java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar SortSam VALIDATION_STRINGENCY=SILENT \
      INPUT=~{aligned_bam} \
      OUTPUT=pipe1 \
      MAX_RECORDS_IN_RAM=~{sort_bam_max_records_in_ram} \
      SORT_ORDER=queryname | \
      java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/picard.jar MergeBamAlignment VALIDATION_STRINGENCY=SILENT \
      ALIGNED_BAM=pipe1 \
      UNMAPPED_BAM=~{unaligned_bam} \
      OUTPUT=pipe2 \
      REFERENCE_SEQUENCE=~{genome_fasta} \
      INCLUDE_SECONDARY_ALIGNMENTS=false \
      PAIRED_RUN=false \
      CLIP_ADAPTERS=false | \
      java -Dsamjdk.compression_level=1 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithInterval VALIDATION_STRINGENCY=SILENT \
      I=pipe2 \
      O=pipe3 \
      INTERVALS=~{gene_intervals} \
      TAG=XG | \
      java -Dsamjdk.compression_level=2 -Xms3000m -jar /software/Drop-seq_tools/jar/dropseq.jar TagReadWithGeneFunction VALIDATION_STRINGENCY=SILENT \
      INPUT=pipe3 \
      O=~{sample_id}_aligned_tagged.bam \
      ANNOTATIONS_FILE=~{refflat}
	"""
def detect_bead_substitution_errors():
	"""
	Identifies cell barcodes with
	aberrant "fixed" UMI bases. If only the last UMI base is fixed as a T, the cell barcode is corrected (the
	last base is trimmed off) and all cell barcodes with identical sequence at the first 11 bases are merged
	together. If any other UMI base is fixed, the reads with that cell barcode are discarded.
	"""
	"""
	java -Dsamjdk.compression_level=1 -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSynthesisErrors VALIDATION_STRINGENCY=SILENT \
      INPUT=~{input_bam} \
      MIN_UMIS_PER_CELL=20 \
      OUTPUT_STATS=~{sample_id}_synthesis_error_stats.txt \
      SUMMARY=~{sample_id}_synthesis_error_summary.txt \
      REPORT=~{sample_id}_synthesis_error_report.txt \
      OUTPUT=temp.bam

    java -Dsamjdk.compression_level=2 -Xmx~{memory} -jar /software/Drop-seq_tools/jar/dropseq.jar DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
      INPUT=temp.bam \
      OUTPUT="~{sample_id}_aligned_tagged_repaired.bam" \
      MIN_UMIS_PER_CELL=20 \
      OUTPUT_REPORT="~{sample_id}_substitution_error_report.txt"
	"""

def bam_tag_histogram(bam: Path, tag: str) -> Path:
	"""
	A key question to answer for your data set is how many cells you want to extract from your BAM. One
	way to estimate this is to extract the number of reads per cell, then plot the cumulative distribution of
	reads and select the “knee” of the distribution.
	We provide a tool to extract the reads per cell barcode in the Drop-seq software called
	BAMTagHistogram. This extracts the number of reads for any BAM tag in a BAM file, and is a general
	purpose tool you can use for a number of purposes. For this purpose, we extract the cell tag “XC”
	"""