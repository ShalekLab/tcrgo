import pandas as pd
import pysam
#import matplotlib
from tcrgo.bam import BAMDict, Read, UMI, Barcode
import tcrgo.io as io
from collections import Counter, defaultdict
from numpy import std

from typing import List, Dict
DataFrame = pd.core.frame.DataFrame
Series = pd.core.series.Series

def subtraction(row: Series, sheet: DataFrame, key: str):
	if row.name not in sheet.index:
		return ''
	variable = sheet.at[row.name, key]
	control = row[key]
	if variable == "None":
		variable = 0
	if control == "None":
		control = 0
	return float(variable) - float(control)

def difference(row: Series, sheet: DataFrame, key: str):
	if row.name not in sheet.index:
		return ''
	print(sheet.at[row.name, key], row[key])
	if sheet.at[row.name, key] == row[key]:
		return "Same"
	else:
		return sheet.at[row.name, key]

def compare(info_control: str, infos_variable: Dict[str, str]):
	control = pd.read_csv(info_control, sep='\t', header=0)
	control.set_index(control["BC"]+control["UMI"])
	for name_sheet, path in infos_variable.items():
		sheet = pd.read_csv(path, sep='\t', header=0)
		sheet.set_index(sheet["BC"]+sheet["UMI"])
		for col in ("nReads", "topVJ_nReads", "topVJ_region", "topVJ_freq", "CDR3_nuc", "CDR3_nReads", "CDR3_freq"):
			name_series = f"{name_sheet}_{col}"
			if col == "topVJ_region" or col == "CDR3_nuc":
				series = control.apply(func=difference, args=(sheet, col), axis=1)
			else:
				series = control.apply(func=subtraction, args=(sheet, col), axis=1)
			control.insert(len(control.columns), name_series, series)
	control.reset_index(drop=True)
	control.to_csv("COMPARISONS.tsv", sep='\t', header=True, index=False)

def id_in_read(bam):
	count = 0
	count_total = 0
	pysam.index(bam)
	bam = pysam.AlignmentFile(bam, "rb")
	for alignment in bam.fetch():
		id = alignment.get_tag("XC") + alignment.get_tag("XU")
		for i in range(0, len(id)-6):
			if alignment.query_sequence.startswith(id[i:]):
				print('\t'.join([alignment.reference_name, alignment.query_name, alignment.query_sequence, id, id[i:], str(i)]))
				count += 1
				break
		count_total += 1
	print(count, count_total, count/count_total)
	"""
	for region in ("TRAC", "TRBC"):
		for alignment in bam.fetch(region):
			id = alignment.get_tag("XC") + alignment.get_tag("XU")
			for i in range(0, len(id)-5):
				if alignment.query_sequence.startswith(id[i:]):
					print('\t'.join([alignment.reference_name, alignment.query_name, alignment.query_sequence, id, id[i:], str(i)]))
	"""

# do you think you could prepare some examples of umis and bc's whose top V/J's are changing 
# and the actual sequences of the reads that are mapping to these different V/J's 
# oh o(so we can make sure it is not just a bad reference genome thing)?
#importlib.reload(ex); ex.analyze_bcumi("/Users/jggatter/Desktop/Projects/tcrgo/out_repaired/queries/190923ShaA_B_raw_repaired_sorted.bam", "bcumis_select.txt", "/Users/jggatter/Desktop/Projects/tcrgo/out_repaired/queries/queries1.txt", "/Users/jggatter/Desktop/Projects/tcrgo/data/macFasCDR3bases_v4.txt", "/Users/jggatter/Desktop/Projects/tcrgo/data/CynoTCRv4.fa")
def analyze_bcumi(bam, ids_select, id_queries, cdr3_file, fasta):
	id_queries = open(id_queries, 'r').read().splitlines()
	ids_select = open(ids_select, 'r').read().splitlines()
	bam = BAMDict(pysam.AlignmentFile(bam, 'rb'))
	pysam.samtools.faidx(str(fasta))
	fasta = pysam.FastaFile(fasta)
	cdr3_positions = io.read_cdr3_file(cdr3_file)

	print("Building BAM Index in memory for fetching VJ queries")
	bam_indexed = pysam.IndexedReads(bam.bam)
	bam_indexed.build()
	print("Built index in memory for fast retrieval")

	for id in ids_select:
		barcode_select = id[:12]
		umi_select = id[12:]
		for id_query in id_queries:
			barcode, umi, query_name = id_query.split('|')
			if barcode == barcode_select and umi == umi_select:
				read = Read(query_name)
				read.parse_alignments(bam_indexed)
				if read.top_V is None or read.top_J is None:
					print(f"Should have gotten a top V and J for {query_name} but did not!")
				for alignment in read:
					read.scores[alignment.reference_name] = alignment.get_tag("AS")
					read.cigars[alignment.reference_name] = alignment.cigarstring
					read.edit_distances[alignment.reference_name] = alignment.get_tag("NM")
					read.phreds[alignment.reference_name] = alignment.get_tag("UQ")
					read.mapqs[alignment.reference_name] = alignment.mapping_quality
				read.get_cdr3_positions(cdr3_positions)
				read.ref_seq_V = fasta.fetch(read.top_V.reference_name)
				read.ref_seq_J = fasta.fetch(read.top_J.reference_name)
				read.get_cdr3_sequence()
				if barcode not in bam:
					bam[barcode] = Barcode(barcode)
				if umi not in bam[barcode]:
					bam[barcode][umi] = UMI(umi)
				bam[barcode][umi][query_name] = read
	
	for barcode_seq, barcode in bam.items():
		for umi_seq, umi in barcode.items():
			umi.find_top_VJ()
			umi.resolve_tcr_identity()
			with open(f"{barcode_seq}_{umi_seq}.tsv", 'w') as umi_file:
				umi_file.write("QNAME\tSEQ\tAlignment Scores\tCIGARs\tEdit Distances\tPhred Scores\ttopV\ttopJ\tisUMItopV\tisUMItopJ\tCDR3\tisUMItopCDR3\n")
				for read in umi.get_reads():
					scores = sorted(read.scores.items(), key=lambda item: item[0], reverse=False)
					cigars = sorted(read.cigars.items(), key=lambda item: item[0], reverse=False)
					edit_distances = sorted(read.edit_distances.items(), key=lambda item: item[0], reverse=False)
					phreds = sorted(read.phreds.items(), key=lambda item: item[0], reverse=False)
					umi_topV, umi_topJ = umi.top_VJ.split('|')
					is_topV = read.top_V.reference_name == umi_topV
					is_topJ = read.top_J.reference_name == umi_topJ
					is_topCDR3 = read.cdr3 == umi.top_cdr3
					umi_file.write(
						f"{read.top_V.query_name}\t{read.top_V.query_sequence}\t"
						f"{scores}\t{cigars}\t{edit_distances}\t{phreds}\t"
						f"{read.top_V.reference_name}\t{read.top_J.reference_name}\t"
						f"{is_topV}\t{is_topJ}\t{read.cdr3}\t{is_topCDR3}\n"
					)

	for barcode_seq, barcode in bam.items():
		for umi_seq, umi in barcode.items():
			scores = defaultdict(list)
			edit_distances = defaultdict(list)
			phreds = defaultdict(list)
			counts = Counter()
			counts_top = Counter()
			with open(f"{barcode_seq}_{umi_seq}.tsv", 'a') as umi_file:
				for read in umi.get_reads():
					print("read.scores", read.scores)
					for ref, score in read.scores.items():
						print(ref, score)
						scores[ref].append(score)
						print(scores[ref])
					for ref, edit_distance in read.edit_distances.items():
						edit_distances[ref].append(edit_distance)
					for ref, phred in read.phreds.items():
						phreds[ref].append(phred)
					for alignment in read:
						counts[alignment.reference_name] += 1
						if not alignment.reference_name in counts_top:
							counts_top[alignment.reference_name] = 0
					counts_top[read.top_V.reference_name] += 1
					counts_top[read.top_J.reference_name] += 1
				references = [ref for ref,count in sorted(counts_top.items(), key=lambda item: item[1], reverse=True)]
				print(f"{references=}")
				print(f"{scores=}")
				print(f"{edit_distances=}")
				print(f"{phreds=}")
				print(f"{counts=}")
				print(f"{counts_top=}")
				umi_file.write("\nV ref\tAlignment Counts\tTop Counts\tAvg. Score\tStdev. Score\tAvg. Edit Distance\tStdev. ED\t Avg. Phred\tStdev. Phred\n")
				for ref in references:
					if not "TRAV" in ref and not "TRBV" in ref:
						continue
					score_avg = 0
					score_stdev = 0
					ed_avg = 0
					ed_stdev = 0
					phred_avg = 0
					phred_stdev = 0
					if scores[ref]:
						score_avg = sum(scores[ref]) / len(scores[ref])
						score_stdev = std(scores[ref])
					if edit_distances[ref]:
						ed_avg = sum(edit_distances[ref]) / len(edit_distances[ref])
						ed_stdev = std(edit_distances[ref])
					if phreds[ref]:
						phred_avg = sum(phreds[ref]) / len(phreds[ref])
						phred_stdev = std(phreds[ref])
					umi_file.write(
						f"{ref}\t{counts[ref]}\t{counts_top[ref]}\t{score_avg:.2f}\t{score_stdev:.2f}\t"
						f"{ed_avg:.2f}\t{ed_stdev:.2f}\t{phred_avg:.2f}\t{phred_stdev:.2f}\n"
					)
				umi_file.write("\nJ ref\tAlignment Counts\tTop Counts\tAvg. Score\tStdev. Score\tAvg. Edit Distance\tStdev. ED\t Avg. Phred\tStdev. Phred\n")
				for ref in references:
					if not "TRAJ" in ref and not "TRBJ" in ref:
						continue
					score_avg = 0
					score_stdev = 0
					ed_avg = 0
					ed_stdev = 0
					phred_avg = 0
					phred_stdev = 0
					if scores[ref]:
						score_avg = sum(scores[ref]) / len(scores[ref])
						score_stdev = std(scores[ref])
					if edit_distances[ref]:
						ed_avg = sum(edit_distances[ref]) / len(edit_distances[ref])
						ed_stdev = std(edit_distances[ref])
					if phreds[ref]:
						phred_avg = sum(phreds[ref]) / len(phreds[ref])
						phred_stdev = std(phreds[ref])
					umi_file.write(
						f"{ref}\t{counts[ref]}\t{counts_top[ref]}\t{score_avg:.2f}\t{score_stdev:.2f}\t"
						f"{ed_avg:.2f}\t{ed_stdev:.2f}\t{phred_avg:.2f}\t{phred_stdev:.2f}\n"
					)
				umi_file.write("\nC ref\tAlignment Counts\tTop Counts\tAvg. Score\tStdev. Score\tAvg. Edit Distance\tStdev. ED\t Avg. Phred\tStdev. Phred\n")
				for ref in references:
					if not "TRAC" in ref and not "TRBC" in ref:
						continue
					score_avg = 0
					score_stdev = 0
					ed_avg = 0
					ed_stdev = 0
					phred_avg = 0
					phred_stdev = 0
					if scores[ref]:
						score_avg = sum(scores[ref]) / len(scores[ref])
						score_stdev = std(scores[ref])
					if edit_distances[ref]:
						ed_avg = sum(edit_distances[ref]) / len(edit_distances[ref])
						ed_stdev = std(edit_distances[ref])
					if phreds[ref]:
						phred_avg = sum(phreds[ref]) / len(phreds[ref])
						phred_stdev = std(phreds[ref])
					umi_file.write(
						f"{ref}\t{counts[ref]}\t{counts_top[ref]}\t{score_avg:.2f}\t{score_stdev:.2f}\t"
						f"{ed_avg:.2f}\t{ed_stdev:.2f}\t{phred_avg:.2f}\t{phred_stdev:.2f}\n"
					)
	for barcode_seq, barcode in bam.items():
		for umi_seq, umi in barcode.items():
			VJ_counts = Counter()
			VJ_cdr3s = defaultdict(set)
			VJ_cdr3_counts = Counter()
			cdr3_counts = Counter()
			with open(f"{barcode_seq}_{umi_seq}.tsv", 'a') as umi_file:
				umi_file.write("\ntopVJ\ttopVJ Count\ttopVJ Freq\tCDR3nuc\tCDR3len\ttopVJ:CDR3 Count\ttopVJ:CDR3 Freq\tAggregated CDR3 Count\tUnaggregated CDR3 Freq\tAggregated CDR3 Freq\n")
				for read in umi.get_reads():
					top_VJ = f"{read.top_V.reference_name}|{read.top_J.reference_name}"
					VJ_counts[top_VJ] += 1
					VJ_cdr3s[top_VJ].add(read.cdr3)
					VJ_cdr3_counts[f"{top_VJ}|{read.cdr3}"] += 1
					cdr3_counts[read.cdr3] += 1
				sum_top_VJs = sum(VJ_counts.values())
				sum_VJ_cdr3s = sum(VJ_cdr3_counts.values())
				sum_cdr3s = sum(cdr3_counts.values())
				VJ_cdr3_sorted = [k.split('|')[-1] for k,v in sorted(VJ_cdr3_counts.items(), key=lambda item: item[1], reverse=True)]
				for top_VJ in [top_VJ for top_VJ,count in sorted(VJ_counts.items(), key=lambda item: item[1], reverse=True)]:
					cdr3s_sorted = [cdr3 for cdr3 in VJ_cdr3_sorted if cdr3 in VJ_cdr3s[top_VJ]]
					for cdr3 in cdr3s_sorted:
						top_VJ_cdr3 = f"{top_VJ}|{cdr3}"
						#sum_VJ_cdr3 = sum(VJ_cdr3_counts[top_VJ_cdr3])
						umi_file.write(
							f"{top_VJ}\t{VJ_counts[top_VJ]}\t{VJ_counts[top_VJ]/sum_top_VJs:.3f}\t{cdr3}\t{len(cdr3)}\t"
							f"{VJ_cdr3_counts[top_VJ_cdr3]}\t{VJ_cdr3_counts[top_VJ_cdr3]/VJ_counts[top_VJ]:.3f}\t"
							f"{VJ_cdr3_counts[top_VJ_cdr3]/sum_VJ_cdr3s:.3f}\t{cdr3_counts[cdr3]}\t{cdr3_counts[cdr3]/sum_cdr3s:.3f}\n"
						)

	print("DONE")

def main():
	infos_variable = {
		"BCip":"/Users/jggatter/Desktop/Projects/tcrgo/out_collapsebcinplace/summary/aggregated_cdr3_info.tsv",
		"BCipU":"/Users/jggatter/Desktop/Projects/tcrgo/out_collapsebcinplaceumi/summary/aggregated_cdr3_info.tsv",
		"BCt":"/Users/jggatter/Desktop/Projects/tcrgo/out_collapsebctag/summary/aggregated_cdr3_info.tsv",
		"BCtU":"/Users/jggatter/Desktop/Projects/tcrgo/out_collapsebctagumi/summary/aggregated_cdr3_info.tsv"
	}
	compare("/Users/jggatter/Desktop/Projects/tcrgo/out_repaired/summary/aggregated_cdr3_info.tsv", infos_variable)

if __name__ == "__main__": main()