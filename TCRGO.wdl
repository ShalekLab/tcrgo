version 1.0

workflow TCRGO {
	input {
		String sample_name
		Array[File] data
		Boolean run_alignment = true
		File fasta
		File cdr3_positions
		Int workers = 32

		String docker = "shaleklab/tcrgo:latest"
		String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		Int preemptible = 2
	}
	if (run_alignment) {
		call preprocessing_and_alignment as alignment {
			input:
				data = data,
				fasta = fasta,
				sample_name = sample_name,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
	}
	call filter_queries {
		input:
			bam = if run_alignment then select_first([alignment.bam_repairedsorted]) else bam,
			workers = workers,
			docker = docker,
			zones = zones,
			preemptible = preemptible
	}
	scatter (query_list in filter_queries.query_list) {
		call reconstruct_tcrs {
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
			cdr3_infos = reconstruct_tcrs.cdr3_info,
			tiebreaks_alignments = select_all(reconstruct_tcrs.tiebreaks_alignments),
			sample_name = sample_name,
			docker = docker,
			zones = zones,
			preemptible = preemptible
	}
	output {
		File aggregated_cdr3_infos = summary.aggregated_cdr3_info
		File aggregated_tiebreaks_alignments = select_all(summary.aggregated_tiebreaks_alignments)
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
		Int task_memory_GB = 32
		String disks = "local-disk 512 HDD"
		Int boot_disk_size_GB = 10
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
		File bam_repairedsorted = "/cromwell_root/out/~{sample_name}_repaired_sorted.bam"
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
		Int number_cpu_threads = 4
		Int task_memory_GB = 16
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 10
	}
	command <<<
		set -e
		cd /scripts/
		python -m filter_queries \
			--minimum_reads ~{minimum_reads} \
			--maximum_reads ~{maximum_reads} \
			--seed ~{seed} \
			--output-path /cromwell_root/out/ \
			--workers ~{workers} \
			~{bam}
	>>>
	output {
		Array[File] query_list = glob("/cromwell_root/out/queries[0-9]*.txt")
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

task reconstruct_tcrs {
	input {
		File bam
		File query_list
		File cdr3_positions
		File fasta
		Float minimum_frequency = 0.3
		Int minimum_cdr3s = 5

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 4
		Int task_memory_GB = 16
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 10
	}
	command <<<
		set -e
		cd /scripts/
		python -m reconstruct_tcrs \
			--fasta ~{fasta} \
			--cdr3-positions-file ~{cdr3_positions} \
			--minimum-frequency ~{minimum_frequency} \
			--minimum-cdr3s ~{minimum_cdr3s} \
			--output-path /cromwell_root/out/ \
			--query-list ~{query_list} \
			~{bam}
	>>>
	output {
		File cdr3_info = glob("/cromwell_root/out/cdr3_info*.tsv")[0]
		File? tiebreaks_alignments = glob("/cromwell_root/out/tiebreaks_alignments*.tsv")[0]
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
		Boolean string_index = false

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 4
		Int task_memory_GB = 16
		String disks = "local-disk 128 HDD"
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
			ALL
		
		# Rename the summary files to contain sample_name
		mv ${output_directory}aggregated_cdr3_info.tsv ${output_directory}~{sample_name}_cdr3_info.tsv
		mv ${output_directory}aggregated_tiebreaks_alignments.tsv ${output_directory}~{sample_name}_tiebreaks_alignments.tsv
	>>>
	output {
		File aggregated_cdr3_info = "/cromwell_root/out/~{sample_name}_cdr3_info.tsv"
		File? aggregated_tiebreaks_alignments = glob("/cromwell_root/out/~{sample_name}_tiebreaks_alignments.tsv")[0]
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