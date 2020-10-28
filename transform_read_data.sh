#!/bin/bash

#The goal of the script is to first transform the Read-1-index-1 format to Read-1-Read-2 format. that way, the TCR CDR3
#sequence_data should be in Read 1, and the associated cell barcode/UMI should be in Read 2.

# First variable is the sequence_data we put in. Second one is the sample name.
# Ex: ./transform_read_data.sh fastq sampleA

sequence_data="$1"
sample="$2"

set -euo pipefail

echo "WARNING: This script requires GNU grep."
echo "On macOS, this can be installed via brew with the following command: brew install grep"
echo "Then any grep commands must be substituted with ggrep."
echo "In future versions of the pipeline this script will be rewritten in Python."

echo Generating ${sample}_IndexReads.txt
awk 'NR%4==2' "${sequence_data}" \
	| rev \
	| tr ATCG TAGC \
	> "${sample}_IndexReads.txt"

echo Generating ${sample}_QSeq.txt
awk 'NR%4==0' "${sequence_data}" > "${sample}_QSeq.txt"

echo Generating ${sample}_BCSeq.txt
awk 'NR%4==1' "${sequence_data}" \
	| grep -o "[ATCGN]*" \
	> "${sample}_BCSeq.txt"

echo Generating ${sample}_QRead2.txt
sed 's~[ATGCN]~@~g' "${sample}_BCSeq.txt" > "${sample}_QRead2.txt"

echo Generating ${sample}_seqHeaders.txt
awk 'NR%4==1' "${sequence_data}" \
	| grep -o "^.*#" \
	> "${sample}_seqHeaders.txt" 
echo Generating ${sample}_seqHeadersRead1.txt
sed 's~#~#/1~' "${sample}_seqHeaders.txt" > "${sample}_seqHeadersRead1.txt"

echo Generating ${sample}_seqHeadersRead2.txt
sed 's~#~#/2~' "${sample}_seqHeaders.txt" > "${sample}_seqHeadersRead2.txt"

echo Generating ${sample}_qualHeadersRead2.txt
sed 's~@~+~' "${sample}_seqHeadersRead2.txt" > "${sample}_qualHeadersRead2.txt"

echo Generating ${sample}_qualHeadersRead1.txt
sed 's~@~+~' "${sample}_seqHeadersRead1.txt" > "${sample}_qualHeadersRead1.txt"

echo Generating ${sample}_TCR.fastq
paste -d '\n' \
	"${sample}_seqHeadersRead2.txt" \
	"${sample}_IndexReads.txt" \
	"${sample}_qualHeadersRead2.txt" \
	"${sample}_QSeq.txt" \
	> "${sample}_TCR.fastq"

echo Generating ${sample}_R1.fastq
paste -d '\n' \
	"${sample}_seqHeadersRead1.txt" \
	"${sample}_BCSeq.txt" \
	"${sample}_qualHeadersRead1.txt" \
	"${sample}_QRead2.txt" \
	> "${sample}_R1.fastq"
