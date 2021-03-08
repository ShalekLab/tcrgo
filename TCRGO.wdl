version 1.0

workflow TCRGO {
	input {
		String sample_name
		Array[File] sequence_data
		
		Boolean run_alignment = true
		File fasta
		File? cdr3_positions
		Int workers = 16

		String docker = "shaleklab/tcrgo:latest"
		String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		Int preemptible = 2
	}
	if (run_alignment) {
		call preprocessing_and_alignment as alignment {
			input:
				data = sequence_data,
				fasta = fasta,
				sample_name = sample_name,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
	}
	call filter_queries {
		input:
			bam = select_first([alignment.bam_repairedsorted, sequence_data[0]]),
			workers = workers,
			docker = docker,
			zones = zones,
			preemptible = preemptible
	}
	scatter (query_list in filter_queries.query_list) {
		call recover_cdr3s {
			input:
				bam = filter_queries.bam_sorted,
				fasta = fasta,
				cdr3_positions = cdr3_positions,
				query_list = query_list,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
	}
	call summary {
		input:
			cdr3_infos = recover_cdr3s.cdr3_info,
			tiebreaks_alignments = recover_cdr3s.tiebreaks_alignments,
			sample_name = sample_name,
			docker = docker,
			zones = zones,
			preemptible = preemptible
	}
	output {
		File? bam_repairedsorted = alignment.bam_repairedsorted
		File aggregated_cdr3_infos = summary.aggregated_cdr3_info
		File aggregated_tiebreaks_alignments = summary.aggregated_tiebreaks_alignments
		File? aggrcollapsed_cdr3_info = summary.aggrcollapsed_cdr3_info
		File? dgraphs = summary.dgraphs
	}
}

task preprocessing_and_alignment {
	input {
		Array[File] data
		File fasta
		String sample_name

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 4
		Int task_memory_GB = 16
		String disks = "local-disk 256 HDD"
		Int boot_disk_size_GB = 15
	}
	command <<<
		set -e
		cd /scripts/
		python -m alignment \
			--fasta ~{fasta} \
			--dropseq /software/dropseq/jar/dropseq.jar \
			--picard /software/dropseq/3rdParty/picard/picard.jar \
			--basename ~{sample_name} \
			--output-path /cromwell_root/out/ \
			~{sep=' ' data}
	>>>
	output {
		File bam_repairedsorted = "/cromwell_root/out/~{sample_name}_repairedsorted.bam"
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}

task filter_queries {
	input {
		File bam
		Int minimum_reads = 5
		Int maximum_reads = 1000
		Int seed = 2020
		Int workers

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 1
		Int task_memory_GB = 8
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 10
	}
	command <<<
		set -e
		cd /scripts/
		python -m filter_queries \
			--minimum-reads ~{minimum_reads} \
			--maximum-reads ~{maximum_reads} \
			--seed ~{seed} \
			--output-path /cromwell_root/out/ \
			--workers ~{workers} \
			~{bam}
	>>>
	output {
		Array[File]+ query_list = glob("/cromwell_root/out/queries[0-9]*.tsv")
		File bam_sorted = bam
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}

task recover_cdr3s {
	input {
		File bam
		File query_list
		File? cdr3_positions
		Boolean is_zero_indexed = false
		Boolean exclude_relatives = false
		File fasta
		Float minimum_frequency = 0.3
		Int minimum_cdr3s = 5

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 1
		Int task_memory_GB = 16
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 10
	}
	command <<<
		set -e
		cd /scripts/
		python -m recover_cdr3s \
			--fasta ~{fasta} \
			~{"--cdr3-positions-file " + cdr3_positions} \
			~{true="--zero-indexed" false='' is_zero_indexed} \
			~{true="--exclude-relatives" false='' exclude_relatives} \
			--minimum-frequency ~{minimum_frequency} \
			--minimum-cdr3s ~{minimum_cdr3s} \
			--output-path /cromwell_root/out/ \
			--query-list ~{query_list} \
			~{bam}
	>>>
	output {
		File cdr3_info = glob("/cromwell_root/out/cdr3_info*.tsv")[0]
		File tiebreaks_alignments = glob("/cromwell_root/out/tiebreaks_alignments*.tsv")[0]
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}

task summary {
	input {
		Array[File] cdr3_infos
		Array[File] tiebreaks_alignments
		String sample_name
		Boolean collapse = true
		Int threshold = 10
		Boolean plot = false
		Boolean string_index = false
		
		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 1
		Int task_memory_GB = 4
		String disks = "local-disk 64 HDD"
		Int boot_disk_size_GB = 10
	}
	command <<<
		set -e
		output_directory="/cromwell_root/out/"
		mkdir -p $output_directory
		cd $output_directory
		
		# Move the files to the PWD
		python <<CODE
		import os
		cdr3_infos = "~{sep=',' cdr3_infos}".split(',')
		tiebreaks_alignments = "~{sep=',' tiebreaks_alignments}".split(',')
		for file in cdr3_infos:
			os.rename(file, os.path.basename(file))
		for file in tiebreaks_alignments:
			os.rename(file, os.path.basename(file))
		CODE
		
		# Run the summary script
		cd /scripts/
		python -m summary \
			~{true="--string-index" false='' string_index} \
			--input-path $output_directory \
			--output-path $output_directory \
			~{true="--collapse" false='' collapse} \
			--threshold ~{threshold} \
			~{true="--plot" false='' plot} \
			ALL
		
		# Rename the summary files to contain sample_name
		mv ${output_directory}aggregated_cdr3_info.tsv \
			${output_directory}~{sample_name}_cdr3_info.tsv
		touch ${output_directory}aggregated_tiebreaks_alignments.tsv
		mv ${output_directory}aggregated_tiebreaks_alignments.tsv \
			${output_directory}~{sample_name}_tiebreaks_alignments.tsv
		touch ${output_directory}aggrcollapsed_cdr3_info.tsv
		mv ${output_directory}aggrcollapsed_cdr3_info.tsv \
			${output_directory}~{sample_name}_aggrcollapsedcdr3_info.tsv
		mkdir -p ${output_directory}dgraphs
		tar -zcf ${output_directory}dgraphs.tar.gz ${output_directory}dgraphs/
	>>>
	output {
		File aggregated_cdr3_info = "/cromwell_root/out/~{sample_name}_cdr3_info.tsv"
		File? aggrcollapsed_cdr3_info = "/cromwell_root/out/~{sample_name}_aggrcollapsedcdr3_info.tsv"
		File? dgraphs = "/cromwell_root/out/dgraphs.tar.gz"
		File aggregated_tiebreaks_alignments = "/cromwell_root/out/~{sample_name}_tiebreaks_alignments.tsv"
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}