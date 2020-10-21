version 1.0

workflow TCRGO {
	input {
		File sample_sheet
		File fasta
		File cdr3_positions
		Int workers = 32

		String docker = "shaleklab/tcrgo:latest"
		String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		Int preemptible = 2
	}
	scatter (sample in read_objects(sample_sheet)) {
		String sample_name = sample.Sample
		String bam_raw = sample.BAM
		call preprocessing_and_alignment as alignment {
			input:
				bam_raw = bam_raw,
				fasta = fasta,
				sample_name = sample_name,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
		call filter_queries {
			input:
				bam = alignment.bam_repaired,
				workers = workers,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
		scatter (file in filter_queries.query_list) {
			call reconstruct_tcrs {
				input:
					bam = alignment.bam_repaired,
					fasta = fasta,
					cdr3_positions = cdr3_positions,
					query_list = file,
					docker = docker,
					zones = zones,
					preemptible = preemptible
			}
		}
		call summary {
			input:
				cdr3_infos = reconstruct_tcrs.cdr3_info,
				sample_name = sample_name,
				docker = docker,
				zones = zones,
				preemptible = preemptible
		}
	}
	output {
		Array[File] aggregated_cdr3_infos = summary.aggregated_cdr3_info
	}
}

task preprocessing_and_alignment {
	input {
		File bam_raw
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
		python -m /tcrgo/alignment \
			--fasta ~{fasta} \
			--dropseq /software/dropseq/jar/dropseq.jar \
			--picard /software/dropseq/3rdparty/picard.jar \
			--basename ~{sample_name} \
			--output-directory /tcrgo/out/alignment/ \
			~{bam_raw}
	>>>
	output {
		File bam_repaired = "/tcrgo/out/alignment/~{sample_name}_repaired.bam"
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
		python -m /scripts/filter_queries \
			--minimum_reads ~{minimum_reads} \
			--maximum_reads ~{maximum_reads} \
			--seed ~{seed} \
			--output-directory /tcrgo/out/queries/ \
			--workers ~{workers} \
			~{bam}
	>>>
	output {
		Array[File] query_list = glob("/tcrgo/out/queries/queries*.txt")
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
		Int minimum_frequency = 0.3
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
		python -m /scripts/reconstruct_tcrs \
			--fasta ~{fasta} \
			--cdr3-positions-file ~{cdr3_positions} \
			--minimum-frequency ~{minimum_frequency} \
			--minimum-cdr3s ~{minimum_cdr3s} \
			--output-directory /tcrgo/out/cdr3/ \
			--query-list ~{query_list} \
			~{bam}
	>>>
	output {
		File cdr3_info = glob("/tcrgo/out/cdr3/cdr3_info*.tsv")[0]
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
		output_directory="/tcrgo/out/cdr3/"
		mkdir -p $output_directory
		cd $output_directory
		python <<CODE
			import os
			# Move the files to the PWD
			cdr3_infos = "~{sep=',' cdr3_infos}".split(',')
			for cdr3_info in cdr3_infos:
				os.rename(cdr3_info, os.path.basename(cdr3_info))
		CODE
		python -m /scripts/summary \
			~{true="--string-index" false='' string_index}
			--input-path $output_directory \
			--output-directory $output_directory \
			ALL
		mv ${output_directory}aggregated_cdr3_info.tsv ${output_directory}~{sample_name}_cdr3_info.tsv
	>>>
	output {
		File aggregated_cdr3_info = "/tcrgo/out/cdr3/~{sample_name}_cdr3_info.tsv"
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