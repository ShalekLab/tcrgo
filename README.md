# Description

TCR Recovery for Seq-Well single-cell RNA-seq data. More to be added here.  

# Requirements

Hardware: macOS (tested), linux (untested) or Windows (untested). Recommended 16GB RAM and at least 50 GB free space on hard drive.  
Software: Python >= 3.8.5, pysam == 0.16.0.1, pandas == 1.1.2, bowtie2 >= 2.4.1, samtools >= 1.3.1, Java 12.0.1 (recommended), Drop-Seq Tools >= 2.3.0  

# Install instructions

It's best you create a new virtual environment. Example using conda on macOS:  
```bash
conda config --add channels conda-forge --add channels bioconda \
&& conda create -y -n tcrgo python=3.8.5 samtools=1.3.1 bowtie2=2.4.1 pysam=0.16.0.1 pandas=1.1.2 \
&& conda activate tcrgo
```
Then `git clone` this repository to your machine.  
For alignment, you must install Java and `git clone` Drop-Seq Tools by the Broad Institute.  

# Run instructions

`cd` to the repository root.  
For help with any of the scripts, use only the `-h` flag for help! Ex: `python -m alignment -h`  
  
## Prepocess and align  
Call alignment.py, a wrapper for Drop-Seq Tools and bowtie2, to tag, trim, align, and repair an unmapped single-end BAM.  
`python -m alignment -d path/to/dropseq.jar -p path/to/picard.jar -f path/to/VJ_reference.fa path/to/unmapped.bam`  
Or first convert paired-end FASTQs to single-end BAM then preprocess and align:  
`python -m alignment -d path/to/dropseq.jar -p path/to/picard.jar -f path/to/VJ_reference.fa path/to/barcode.fastq.gz path/to/biological.fastq.gz`  
 
## Find and distribute VJ-mapped queries
Call filter_queries.py to write a list of BAM queries (reads) which contain alignments to V and J subregions.  
`python -m filter_queries -lt 5 -gt 1000 /path/to/alignment.bam`  
where alignment.bam is a BAM produced by the alignment workflow.  
The `-w` flag can be used to divide the list of queries up into *w* files for processing by workers on an HPC or the cloud.  
  
## Determine consensus VJ pairs and CDR3 sequences per UMI
Call reconstruct_tcrs.opy to find the top-scoring V and J alignments per read, the VJ pair alignment consensus per UMI, and the consensus CDR3 sequence per UMI.  
`python -m reconstruct_tcrs -f path/to/VJ_reference.fasta -c path/to/VJ_cdr3_positions.tsv -w 1 /path/to/alignment.bam`  
where alignment.bam is a BAM produced by the alignment workflow.  
The `-w` flag can be used to run the script on the *w*th query list produced by filter_queries.py.  
  
## Aggregate statistics from reconstruct_tcrs
Call summary.py to combine all cdr3_info***w***.tsv files produced by reconstruct_tcrs.py into one large tsv file.  
`python -m summary -i path/to/cdr3_info_folder/ 1:w`  
where `w` is the number of the final cdr3_info***w***.tsv file you wish to aggregate.  



