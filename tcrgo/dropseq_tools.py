"""
Subprocess calls to Drop-seq Tools
"""
import subprocess as sp
import pysam
import os
from textwrap import dedent
import math
from typing import List, Tuple, Dict
from tcrgo.log import Log

log = Log("root")

def execute(command: str) -> str:
	try:
		output = sp.check_output(command.split())
	except sp.CalledProcessError:
		log.error(
			f"{command.split()} terminated with non-zero exit status, "
			"please check the above messages for more detail."
		)
	return output

def test_tools(dropseq_jar: str, picard_jar: str, aligner: str):
	"""Notify user if tools could not be called from commandline"""
	if not os.path.isfile(dropseq_jar):
		log.error(f"Could not find {dropseq_jar}")
	if not os.path.isfile(picard_jar):
		log.error(f"Could not find {picard_jar}")
	execute("samtools --version")
	if aligner == "bowtie2":
		execute("bowtie2 --version")
	elif aligner == "star":
		execute("star --version")
	elif aligner == "bwa":
		log.error("BWA is currently unsupported.")

def bam2fq(bam: str, fastq: str) -> str:
	command = f"samtools bam2fq -0 {fastq} {bam}"
	execute(command)
	return fastq

def transform_read_data(fastq_singleend: str, sample_name: str, fastq_barcode: str, fastq_biological: str, output_path: str) -> Tuple[str, str]:
	execute("chmod +x transform_read_data.sh")
	command = f"./transform_read_data.sh {fastq_singleend} {sample_name} {output_path}"
	execute(command)
	return fastq_barcode, fastq_biological

def fastq_to_bam(picard_jar: str, fastq_barcode: str, fastq_biological: str, bam_unmapped: str, sample_name: str) -> str:
	"""Convert FASTQ R1 and R2 to unmapped BAM using Picard's FastqToSam"""
	command = \
		f"""
		java -Dsamjdk.compression_level=6 -Xmx4g -jar {picard_jar} FastqToSam \
			SAMPLE_NAME={sample_name} \
			OUTPUT={bam_unmapped} \
			FASTQ={fastq_barcode} \
			FASTQ2={fastq_biological} \
			QUALITY_FORMAT=Standard \
			SORT_ORDER=queryname \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
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
		java -Xmx4000m -jar {dropseq_jar} \
		TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_untagged} \
			OUTPUT={bam_celltagged} \
			BARCODE_QUALITY_TAG=CY \
			SUMMARY={summary_tagged_cellular} \
			BASE_RANGE=1-12 \
			BASE_QUALITY=10 \
			BARCODED_READ=1 \
			DISCARD_READ=false \
			TAG_NAME=CR \
			NUM_BASES_BELOW_QUALITY=1 \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	log.sep()
	execute(command)
	log.sep()
	command = \
		f"""
		java -Xmx4000m -jar {dropseq_jar} \
			TagBamWithReadSequenceExtended VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_celltagged} \
			OUTPUT={bam_idtagged} \
			BARCODE_QUALITY_TAG=QX \
			SUMMARY={summary_tagged_molecular} \
			BASE_RANGE=13-20 \
			BASE_QUALITY=10 \
			BARCODED_READ=1 \
			DISCARD_READ=true \
			TAG_NAME=RX \
			NUM_BASES_BELOW_QUALITY=1 \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	log.sep()
	execute(command)
	log.sep()
	return bam_idtagged

def filter_bam(dropseq_jar: str, bam_tagged: str, bam_filtered: str):
	execute(
		f"""
		java -Xms4g -jar {dropseq_jar} \
		FilterBam \
			INPUT={bam_tagged} \
			OUTPUT={bam_filtered} \
			TAG_REJECT=XQ \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	)
	# FILTER_PCR_DUPES=true ? what are pcr dupes anyways versus reads of same UMI?
	# SUM_MATCHING_BASES=int 
	# 	Retain reads that have at least this many M bases total in the cigar string. 
	#	This sums all the M's in the cigar string.  Default value: null. 

	return bam_filtered

def trim_bam(dropseq_jar: str, bam_idtagged: str, bam_trimmed: str) -> str:
	dirname = os.path.dirname(bam_trimmed)
	bam_adaptertrimmed = os.path.basename(bam_trimmed).replace("_trimmed.bam", "_adaptertrimmed.bam")
	bam_adaptertrimmed = os.path.join(dirname, bam_adaptertrimmed)
	summary_adaptertrimmed = os.path.join(dirname, "summary_trimming_adapter.txt")
	summary_polyAtrimmed = os.path.join(dirname, "summary_trimming_polyA.txt")
	log.sep()
	execute(
		f"""
		java -Xmx4g -jar {dropseq_jar} 
		TrimStartingSequence VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_idtagged} \
			OUTPUT={bam_adaptertrimmed} \
			OUTPUT_SUMMARY={summary_adaptertrimmed} \
			SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
			MISMATCHES=0 \
			NUM_BASES=5 \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	)
	log.sep()
	execute(
		f"""
		java -Xmx4g -jar {dropseq_jar}
		PolyATrimmer VALIDATION_STRINGENCY=SILENT \
			INPUT={bam_adaptertrimmed} \
			OUTPUT={bam_trimmed} \
			OUTPUT_SUMMARY={summary_polyAtrimmed} \
			ADAPTER=^RX^CRACGTACTCTGCGTTGCTACCACTG \
			MISMATCHES=0 \
			NUM_BASES=6 \
			USE_NEW_TRIMMER=true \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	)
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
	

def bowtie2(bam_trimmed: str, fasta: str, sam_aligned: str) -> str:
	"""Sort BAM using Picard then align to reference using bowtie2"""
	if not sam_aligned.endswith("_aligned.sam"):
		log.error("sam_aligned must end in '_aligned.sam'.")
	bam_sorted = bam_trimmed.replace("_trimmed.bam", "_sorted.bam")
	pysam.sort("-n", "-o", bam_sorted, bam_trimmed) # Sort by QNAME
	pysam.index(bam_sorted)

	index_prefix = os.path.basename(fasta.split('.')[0])
	index_prefix = os.path.join(os.path.dirname(sam_aligned), index_prefix)
	command = f"bowtie2-build -f {fasta} --threads 8 {index_prefix}"
	log.sep()
	execute(command)
	command = \
		f"""
		bowtie2 -a --very-sensitive-local --norc \
			-x {index_prefix} \
			-p 8 \
			-b {bam_sorted} \
			-S {sam_aligned} \
		"""
	execute(command)
	log.sep()
	return sam_aligned

def sam_to_fastq(picard_jar: str, bam: str, fastq: str) -> str:
	command = \
		f"""
		java -Xms4g -jar {picard_jar} SamToFastq \
			I={bam}
			FASTQ={fastq}
		"""
	execute(command)
	return fastq

def bwa(fastq_trimmed: str, fasta: str, sam_aligned: str) -> str:
	command = f"bwa index {fasta}"
	execute(command)
	command = f"bwa mem -M -t 8 -o {sam_aligned} -a {fasta} {fastq_trimmed}"
	execute(command)
	return sam_aligned

def star(fastq_trimmed: str, fasta: str, basename: str, sam_aligned: str) -> str:
	path = os.path.dirname(fasta)
	ref_path = os.path.join(path, "ref")
	ref_exists = True
	if not os.path.exists(ref_path):
		os.makedirs(ref_path)
		ref_exists = False
	else:
		log.warn("STAR ref files detected, will not be generating a new reference.")

	genome_length = 0
	with open(fasta, 'r') as fastafile:
		for line in fastafile:
			if not line.startswith('>'):
				genome_length += len(line.strip())
	genomeSAindexNbases = math.floor(min(14, math.log(genome_length, 2) / 2 - 1))

	if not ref_exists:
		log.sep()
		execute(
			f"""
			star --runMode genomeGenerate \
				--runThreadN 4 \
				--genomeDir {ref_path} \
				--genomeFastaFiles {fasta} \
				--genomeSAindexNbases {genomeSAindexNbases}
			"""
		)
	log.sep()
	execute(
		f"""
		star --runMode alignReads \
			--runThreadN 4 \
			--genomeDir {ref_path} \
			--sjdbOverhang 149 \
			--outFilterScoreMinOverLread 0 \
			--outFilterMatchNminOverLread 0 \
			--outFilterMatchNmin 0 \
			--outFilterMultimapNmax 999 \
			--readFilesIn {fastq_trimmed} \
			--outFileNamePrefix {basename}
		"""
	)
	log.sep()
	os.rename(basename+"Aligned.out.sam", sam_aligned)
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
		log.warn(f"{fasta_dict} already exists! Removing and creating new dict.")
		os.remove(fasta_dict)
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
			INCLUDE_SECONDARY_ALIGNMENTS=true \
			MAX_INSERTIONS_OR_DELETIONS=-1 \
			PAIRED_RUN=false \
			SORT_ORDER=coordinate \
			CLIP_ADAPTERS=false \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		""" # There was a q after the last 'false'?
	log.sep()
	execute(command)
	log.sep()
	return merged_bam

def fauxtag_exons(bam_merged: str, bam_exontagged: str) -> str:
	"""
	This isn't really 'faux-tagging' anymore. This function writes a new BAM where
	mapped forward reads are tagged accordingly with either TRA or TRB so that
	Drop-Seq Tools merge and repair steps can proceed and potentially use this
	info for Barcode collapsing.
	"""
	log.sep()
	pysam.index(bam_merged)
	log.info(f"Faux-tagging {bam_merged} to produce {bam_exontagged}.")
	with pysam.AlignmentFile(bam_merged, 'rb') as untagged:
		with pysam.AlignmentFile(bam_exontagged, 'wb', template=untagged) as tagged:
			exon_tags = [("gf", "CODING", "Z"), ("gn", "UKN", "Z"), ("gs", "+", "Z")]
			for query in untagged:
				if query.is_unmapped or query.is_reverse:
					continue
				if "TRA" in query.reference_name:
					exon_tags[1] = ("gn", "TRA", "Z")
				elif "TRB" in query.reference_name:
					exon_tags[1] = ("gn", "TRB", "Z")
				else:
					raise Exception(f"Unknown gene in {query.reference_name}!")
				tags = list(query.get_tags())
				query.set_tags(tags + exon_tags)
				tagged.write(query)
	log.sep()
	return bam_exontagged

def bead_substitution_errors(dropseq_jar: str, bam: str, bam_out: str, min_umis_per_cell: int=1) -> str:
	dirname = os.path.dirname(bam_out)
	report_substitution_error = os.path.join(dirname, "report_substitution_error.txt")
	log.sep()
	execute(
		f"""
		java -Xmx4g -jar {dropseq_jar} \
		DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
			VERBOSITY=DEBUG \
			INPUT={bam} \
			OUTPUT={bam_out} \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			READ_MQ=10 \
			FREQ_COMMON_SUBSTITUTION=0.8 \
			CELL_BARCODE_TAG=CR \
			OUT_CELL_BARCODE_TAG=CR \
			MOLECULAR_BARCODE_TAG=RX \
			NUM_THREADS=4 \
			OUTPUT_REPORT={report_substitution_error} \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	)
	log.sep()
	return

def bead_synthesis_errors(dropseq_jar: str, bam: str, bam_out: str, min_umis_per_cell: int=1) -> str:
	dirname = os.path.dirname(bam_out)
	stats_synthesis_error = os.path.join(dirname, "stats_synthesis_error.txt")
	summary_synthesis_error = os.path.join(dirname, "summary_synthesis_error.txt")
	report_synthesis_error = os.path.join(dirname, "report_synthesis_error.txt")
	log.sep()
	execute(
		f"""
		java -Xmx4g -jar {dropseq_jar} \
		DetectBeadSynthesisErrors VALIDATION_STRINGENCY=SILENT \
			VERBOSITY=DEBUG \
			INPUT={bam} \
			OUTPUT_STATS={stats_synthesis_error} \
			SUMMARY={summary_synthesis_error} \
			REPORT={report_synthesis_error} \
			CELL_BARCODE_TAG=CR \
			MOLECULAR_BARCODE_TAG=RX \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			READ_MQ=10 \
			NUM_THREADS=4 \
			OUTPUT={bam_out} \
			USE_JDK_DEFLATER=true \
			USE_JDK_INFLATER=true
		"""
	)
	log.sep()
	return bam_out

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
	stats_synthesis_error = os.path.join(dirname, "stats_synthesis_error.txt")
	summary_synthesis_error = os.path.join(dirname, "summary_synthesis_error.txt")
	report_substitution_error = os.path.join(dirname, "report_substitution_error.txt")
	report_synthesis_error = os.path.join(dirname, "report_synthesis_error.txt")
	command = \
		f"""
		java -Dsamjdk.compression_level=2 -Xmx4g -jar {dropseq_jar} \
		DetectBeadSubstitutionErrors VALIDATION_STRINGENCY=SILENT \
			VERBOSITY=DEBUG \
			INPUT={bam_exontagged} \
			OUTPUT={bam_subrepaired} \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			READ_MQ=10 \
			CELL_BARCODE_TAG=CR \
			MOLECULAR_BARCODE_TAG=RX \
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
			OUTPUT_STATS={stats_synthesis_error} \
			SUMMARY={summary_synthesis_error} \
			REPORT={report_synthesis_error} \
			CELL_BARCODE_TAG=CR \
			MOLECULAR_BARCODE_TAG=RX \
			MIN_UMIS_PER_CELL={min_umis_per_cell} \
			READ_MQ=10 \
			NUM_THREADS=4 \
			OUTPUT={bam_repaired}
		"""
	log.sep()
	execute(command)
	log.sep()
	bam_repairedsorted = bam_repaired.replace(".bam", "sorted.bam") 
	print(f"Sorting bam to produce {bam_repairedsorted}...")
	pysam.sort("-o", bam_repairedsorted, bam_repaired)
	return bam_repairedsorted

# TODO: Python not integrated into script.
def collapsebarcodesinplace(dropseq_jar: str, bam: str, bam_output: str):
	command = \
		f"""
		java -jar {dropseq_jar} CollapseBarcodesInPlace \
			INPUT={bam} \
			OUTPUT={bam_output} \
			PRIMARY_BARCODE=CR \
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
def collapseumiswithcontext(dropseq_jar: str, bam: str, bam_output: str):
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

# TODO: Python untested.
def collapsebarcodeswithcontext(dropseq_jar: str, bam: str, bam_output: str):
	command = \
		f"""
		java -Xmx16g -jar {dropseq_jar} CollapseTagWithContext \
			VERBOSITY=DEBUG \
			INPUT={bam} \
			COLLAPSE_TAG=XC \
			CONTEXT_TAGS=gf \
			OUT_TAG=XC \
			OUTPUT={bam_output} \
			EDIT_DISTANCE=1 \
			FIND_INDELS=true \
			READ_MQ=10
			MIN_COUNT=1 \
			DROP_SMALL_COUNTS=false \
			NUM_THREADS=12 \
			USE_JDK_INFLATER=true \
			USE_JDK_DEFLATER=true \
			LOW_MEMORY_MODE=true
		""" # COUNT_TAGS option?
	log.sep()
	execute(command)
	log.sep()
	return bam_output
