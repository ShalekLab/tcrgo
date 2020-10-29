#!/bin/bash

#The goal of the script is to first transform the Read-1-index-1 format to Read-1-Read-2 format. that way, the TCR CDR3
#sequence_data should be in Read 1, and the associated cell barcode/UMI should be in Read 2.

# First variable is the sequence_data we put in. Second one is the sample name.
# Ex: ./transform_read_data.sh fastq sampleA

sequence_data="$1"
sample="$2"
output_path="$3"

set -euo pipefail

echo "WARNING: This script requires GNU grep."
echo "On macOS, this can be installed via brew with the following command: brew install grep"
echo "Then any grep commands must be substituted with ggrep."
echo "In future versions of the pipeline this script will be rewritten in Python."

echo Generating ${output_path}${sample}_IndexReads.txt
awk 'NR%4==2' "${sequence_data}" \
	| rev \
	| tr ATCG TAGC \
	> "${output_path}${sample}_IndexReads.txt"

echo Generating ${output_path}${sample}_QSeq.txt
awk 'NR%4==0' "${sequence_data}" > "${output_path}${sample}_QSeq.txt"

echo Generating ${output_path}${sample}_BCSeq.txt
awk 'NR%4==1' "${sequence_data}" \
	| grep -o "[ATCGN]*" \
	> "${output_path}${sample}_BCSeq.txt"

echo Generating ${output_path}${sample}_QRead2.txt
sed 's~[ATGCN]~@~g' "${output_path}${sample}_BCSeq.txt" > "${output_path}${sample}_QRead2.txt"

echo Generating ${output_path}${sample}_seqHeaders.txt
awk 'NR%4==1' "${sequence_data}" \
	| grep -o "^.*#" \
	> "${output_path}${sample}_seqHeaders.txt" 
echo Generating ${output_path}${sample}_seqHeadersRead1.txt
sed 's~#~#/1~' "${output_path}${sample}_seqHeaders.txt" > "${output_path}${sample}_seqHeadersRead1.txt"

echo Generating ${output_path}${sample}_seqHeadersRead2.txt
sed 's~#~#/2~' "${output_path}${sample}_seqHeaders.txt" > "${output_path}${sample}_seqHeadersRead2.txt"

echo Generating ${output_path}${sample}_qualHeadersRead2.txt
sed 's~@~+~' "${output_path}${sample}_seqHeadersRead2.txt" > "${output_path}${sample}_qualHeadersRead2.txt"

echo Generating ${output_path}${sample}_qualHeadersRead1.txt
sed 's~@~+~' "${output_path}${sample}_seqHeadersRead1.txt" > "${output_path}${sample}_qualHeadersRead1.txt"

echo Generating ${output_path}${sample}_TCR.fastq
paste -d '\n' \
	"${output_path}${sample}_seqHeadersRead2.txt" \
	"${output_path}${sample}_IndexReads.txt" \
	"${output_path}${sample}_qualHeadersRead2.txt" \
	"${output_path}${sample}_QSeq.txt" \
	> "${output_path}${sample}_TCR.fastq"

echo Generating ${output_path}${sample}_R1.fastq
paste -d '\n' \
	"${output_path}${sample}_seqHeadersRead1.txt" \
	"${output_path}${sample}_BCSeq.txt" \
	"${output_path}${sample}_qualHeadersRead1.txt" \
	"${output_path}${sample}_QRead2.txt" \
	> "${output_path}${sample}_R1.fastq"
