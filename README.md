# Description

The computational Seq-Well TCR alignment pipeline identifies and quantifies the TCR regions associated with each cell barcode in the sequencing data resulting from TCR enrichment protocol described in Tu et al. 2019. We have reengineered the Seq-Well TCR recovery pipeline originally presented by Tu et al. to make significant enhancements to the user and developer experience while delivering comparable results. Chiefly, we have improved the performance, reporting, tunability, portability, and maintainability of the computational pipeline in a new Python program, TCRGO. This pipeline consists of four main steps, preprocessing and alignment, filtering using alignment information, CDR3 recovery, and summarizing of the data. The preprocessing and alignment step consists of conversion to single-end BAM format, barcode tagging, read trimming, multimap alignment to TRA and TRB gene segments, and barcode error correction. Multimap alignment to V, J, and C TRA and TRB gene segments is performed via Bowtie2 while all other preprocessing tasks are now handled by trusted tools such as Samtools, Picard, and Picard extension Drop-Seq Tools. The filtering step iterates over the alignment data and distributes reads that multimap to both a V and J segment to lists which can be used in parallel by the CDR3 recovery step. The CDR3 recovery step quantifies the alignment identity data for each unique transcript and determines the unique transcript’s most likely CDR3 sequence using a greedy, optimized, heuristic-driven algorithm. The alignment identity and CDR3 information is gathered and tabulated during the summary step. The filtering, recovery, and analysis steps, originally scripted in Bourne Again Shell and Matlab, have been refactored in Python to afford better maintainability of the codebase, reporting for easier troubleshooting, and boosted performance via harnessing highly efficient and powerful packages such as pysam, biopython, numpy, and pandas. Furthermore, the pipeline has been fitted for cloud computing through wrapping it in Docker and the Workflow Description Language. Cloud deployability will increase accessibility to scientists on cloud research platforms such as Terra and enable parallelization of the pipeline’s data processing jobs.

# Requirements

Hardware: macOS (tested), linux (untested), or Windows (untested). Recommended 16GB RAM and at least 64 GB free space on hard drive.  
Software: Python ≈ 3.8.5, pysam ≈ 0.16.0.1, pandas ≈ 1.1.2, bowtie2 >= 2.4.1, samtools ≈ 1.3.1 python-igraph ≈ 0.8.3, biopython = 1.77, Java 12.0.1 (recommended), Drop-Seq Tools >= 2.4.0  

# Instructions for running on the cloud

This program has been wrapped as a Workflow Description Language script, TCRGO.wdl. The latest version of this workflow can be found here on [Dockstore](https://dockstore.org/workflows/github.com/ShalekLab/tcrgo/TCRGO), where it can be exported to an analysis platform of your choice. The Docker image itself is hosted [here on Docker Hub](https://hub.docker.com/r/shaleklab/tcrgo).  

On Terra, scattering of samples is possible through using the [Terra Data Table](https://support.terra.bio/hc/en-us/articles/360025758392-Managing-data-with-workspace-tables-). 

# Instructions for a local installation

Please see the Hardware requirements before considering a local installation.

## Installing

It's strong recommended you create a new virtual environment. Example using `conda` on a Unix-like system:  
```bash
conda config --add channels conda-forge --add channels bioconda \
&& conda create -y -n tcrgo python=3.8.5 samtools=1.3.1 bowtie2=2.4.1 pysam=0.16.0.1 pandas=1.1.2 python-igraph=0.8.3 biopython=1.77 \
&& conda activate tcrgo
```
Then `git clone` this repository to your machine.  
For alignment, you must install Java and `git clone` [Drop-Seq Tools](https://github.com/broadinstitute/Drop-seq) by the Broad Institute. The Picard and Java jar files must be entered as arguments to the alignment script.  

## Running

Change the present working directory (`cd` on Unix-like systems) to the repository root folder.  
For help with any of the scripts, use only the `-h` flag for help! Ex: `python -m alignment -h`  
On a Unix-like system, if there are spaces in any of your pathnames please escape them using `\`! Example: `/path/to/my data/my sample.bam` becomes `/paht/to/my\ data/my\ sample.bam`.
  
### Prepocess and align  
Call alignment.py, a wrapper for Drop-Seq Tools and bowtie2, to tag, trim, align, and barcode-repair a BAM, single-end, or .  
`python -m alignment -d path/to/dropseq.jar -p path/to/picard.jar -f path/to/VJ_reference.fasta -b mysamplename path/to/raw.bam`  
Or first convert paired-end FASTQs to single-end BAM then preprocess and align:  
`python -m alignment -d path/to/dropseq.jar -p path/to/picard.jar -f path/to/VJ_reference.fasta -b mysamplename path/to/barcode.fastq.gz,path/to/biological.fastq.gz`  

Here is an overview of the preprocessing/alignment pipeline:
```
Note: `<basename>` is inputted using the `--basename` argument.
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
```

### Find and distribute VJ-mapped queries
Call filter_queries.py to write a list of BAM queries (reads) which contain alignments to V and J subregions.  
`python -m filter_queries -lt 5 -gt 1000 /path/to/alignment.bam`  
where alignment.bam is a BAM produced by the alignment workflow.  
The `-w` flag can be used to divide the list of queries up into *w* files for processing by workers on an HPC or the cloud.  
  
### Determine consensus VJ pairs and CDR3 sequences per UMI
Call recover_cdr3s.py to find the top-scoring V and J alignments per read, the VJ pair alignment consensus per UMI, and the consensus CDR3 sequence per UMI.  
`python -m recover_cdr3s -f path/to/VJ_reference.fasta -mf 0.3 -mc 5 -q out/queries1.tsv /path/to/alignment.bam`  
where alignment.bam is a BAM produced by the alignment workflow.  

The `--cdr3-positions-file` argument can be used to input a 2-column TSV file containing the names of FASTA TR(A|B)(V|J) segments and that entry's corresponding CDR3 start or ending codon position in the nucleotide sequence. Each entry name in the CDR3 positions file must exactly match its respective entry in the FASTA. If the argument is not specified, the program will attempt to identify the most likely CDR3 start and end positions.

### Aggregate statistics from recover_cdr3s
Call summary.py to combine all cdr3_info***w***.tsv files produced by recover_cdr3s.py into one large tsv file. The `--collapse` mode will attempt to collapse UMIs by their cell and molecular barcodes; an intensive process which consumes a considerable amount of time. 
`python -m summary -i path/to/cdr3_info_folder/ --collapse ALL`



