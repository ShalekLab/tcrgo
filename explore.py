"""This is a notebook file I worked on in Spyder IDE"""

import pandas as pd
import matplotlib.pyplot as plt
import pysam

import sys; sys.path.insert(0, "/Users/jggatter/Desktop/Projects/tcrgo/")
import tcrgo
from tcrgo.bam import Read
import tcrgo.io

# %% Load Data
split_int = lambda s: [int(i) for i in s.split(',')]
split_str = lambda s: s.split(',')
converters={
	"FLAGs":split_int,
	"Regions":split_str,
	"Scores":split_int,
	"CIGARs":split_str,
	"MAPQs":split_int,
	"Phreds":split_int,
	"EditDistance":split_int
}
reads = pd.read_csv(
	"/Users/jggatter/Desktop/Projects/tcrgo/reads.tsv",
	sep='\t',
	header=0,
	index_col="QNAME",
	converters=converters,
)
head = reads.head()
print("LOADED")

# %% Basic charts

def bar(series, title="Title", xlabel="X", ylabel="Y", color='b'):
	plt.figure(figsize=(10, 10))
	counts = series.value_counts(dropna=False)
	plt.bar(counts.index.to_list(), counts, color=color)
	plt.xlabel(xlabel)
	plt.xticks(rotation=90)
	plt.ylabel(ylabel)
	plt.show()
	

def staggered_explode(values, offset=0.2, threshold=0.05):
	explode = []
	offset=0.2
	initial_offset=offset
	for value in values:
		if value <= threshold:
			explode.append(offset)
			if offset == initial_offset:
				offset += offset
			else:
				offset /= 2
		else:
			explode.append(0)
	return explode

def ascending_explode(values, initial_offset=0.05, increment_offset=0.05, threshold=0.05):
	explode = []
	offset = initial_offset
	for value in values:
		if value <= threshold:
			explode.append(offset)
			offset += increment_offset
		else:
			explode.append(0)
	return explode

def standard_explode(values, offset=0.2, threshold=0.05):
	explode = []
	for value in values:
		if value <= threshold:
			explode.append(offset)
		else:
			explode.append(0)
	return explode

def pie(series):
	plt.figure(figsize=(20, 20))
	percentages = series.value_counts(normalize=True, dropna=False, sort=False)
	labels = percentages.index.to_list()
	#explode = staggered_explode(percentages)
	#explode = ascending_explode(percentages)
	explode = standard_explode(percentages)
	fig, ax = plt.subplots()
	ax.pie(percentages, explode=explode, labels=labels, autopct='%1.1f%%', labeldistance=1.5, rotatelabels=True)
	ax.axis('equal')
	plt.tight_layout()
	plt.show()
	
def hist(series, xlabel="X", ylabel="Y", title="Histogram", bins=10, log=False):
	plt.figure(figsize=(10, 10))
	counts = series.value_counts(dropna=False)
	plt.hist(counts, bins, density=True, log=log, facecolor='g', alpha=0.75)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.grid(True)
	plt.show()

# %% Basic

pie(reads["hasVJ"])
pie(reads["CDR3status"])

# %% Flags setup
flags = reads["FLAGs"].explode()
# All these flags are unpaired mapped
# There were no unmapped.
flag_alias = {
	0:"Primary forward strand",
	16:"Primary reverse strand",
	256:"Non-primary forward strand",
	272:"Non-primary reverse strand"
}
flags_aliased = flags.map(flag_alias)

# %% Flags plot
bar(flags_aliased, xlabel="FLAG", ylabel="Count")
pie(flags_aliased)

# %% CDR3 status among forward strands
#flags = reads["FLAGs"].explode()
forward_primary_index = flags.index[flags == 0]
forward_nonprimary_index = flags.index[flags == 256]
forward_index = flags.index[(flags == 0) | (flags == 256)]

forward = reads.loc[forward_index.to_list()]
forward_primary = reads.loc[forward_primary_index.to_list()]
forward_nonprimary = reads.loc[forward_nonprimary_index.to_list()]

# %% Forward CDR3 Plot
pie(reads["CDR3status"])
pie(forward_primary["CDR3status"])
pie(forward_nonprimary["CDR3status"])
pie(forward["CDR3status"])
# For forward mapped reads, approximately 1/3rd had error while constructing CDR3
# Perhaps the other two thirds didn't error by happenstance?

# %% UMIs per barcode, reads per UMI
#For mapped:
# UMIs per barcode histogram
# Reads per UMIs histogram

barcode_umi_counts = reads.groupby(['Barcode', 'UMI']).size().reset_index().rename(columns={0:'count'})
hist(reads["Barcode"], bins=20, log=True, title="non-unique UMIs per Barcode", xlabel="# UMIs", ylabel="Log-scale Proportion of Barcodes")
hist(barcode_umi_counts["Barcode"], bins=20, log=True, title="unique UMIs per Barcode", xlabel="# UMIs", ylabel="Log-scale Proportion of Barcodes")
hist(reads["UMI"], bins=20, log=True, title="Reads per UMIs", xlabel="# Reads", ylabel="Log-scale Proportion of UMIs")

# %%
# For successful CDR3s I wanna see what the Alignment Scores, Phreds, etc. look like versus
# the ERR ones

cdr3_error = reads.loc[reads['CDR3status'] == 'ERR']
cdr3_success = reads.loc[(reads['CDR3status'] == 'COM') | (reads['CDR3status'] == 'INC')]
cdr3_success = cdr3_success.sample(n=1000, replace=False, random_state=999)

# %% 
def get_from_lists(row, subregion_col):
	i = row[subregion_col]
	for col in ("Regions", "Scores", "Phreds", "MAPQs", "EditDistance"):
		row[subregion_col+'_'+col] = row[col][i]
	return row
# Quick
cdr3_error.apply(func=get_from_lists, axis=1, args=("topV",))
cdr3_error.apply(func=get_from_lists, axis=1, args=("topJ",))
head_err = cdr3_error.head()
# %% Takes forever
cdr3_success.apply(func=get_from_lists, axis=1, args=("topV",))
cdr3_success.apply(func=get_from_lists, axis=1, args=("topJ",))
head_suc = cdr3_success.head()

# %%

def get_from_list(row, value, key):
	return row[value][row[key]]

cdr3_error["topV_Region"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Regions", "topV"))
cdr3_error["topJ_Region"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Regions", "topJ"))
cdr3_error["topV_Score"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Scores", "topV"))
cdr3_error["topJ_Score"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Scores", "topJ"))
cdr3_error["topV_Phred"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Phreds", "topV"))
cdr3_error["topJ_Phred"] = cdr3_error.apply(func=get_from_list, axis=1, args=("Phreds", "topJ"))
cdr3_error["topV_MAPQ"] = cdr3_error.apply(func=get_from_list, axis=1, args=("MAPQs", "topV"))
cdr3_error["topJ_MAPQ"] = cdr3_error.apply(func=get_from_list, axis=1, args=("MAPQs", "topJ"))
cdr3_error["topV_EditDis"] = cdr3_error.apply(func=get_from_list, axis=1, args=("EditDistance", "topV"))
cdr3_error["topJ_EditDis"] = cdr3_error.apply(func=get_from_list, axis=1, args=("EditDistance", "topJ"))
head_err = cdr3_error.head()

# %%
def get_from_list(row, value, key):
	return row[value][row[key]]

cdr3_success["topV_Region"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Regions", "topV"))
cdr3_success["topJ_Region"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Regions", "topJ"))
cdr3_success["topV_Score"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Scores", "topV"))
cdr3_success["topJ_Score"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Scores", "topJ"))
cdr3_success["topV_Phred"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Phreds", "topV"))
cdr3_success["topJ_Phred"] = cdr3_success.apply(func=get_from_list, axis=1, args=("Phreds", "topJ"))
cdr3_success["topV_MAPQ"] = cdr3_success.apply(func=get_from_list, axis=1, args=("MAPQs", "topV"))
cdr3_success["topJ_MAPQ"] = cdr3_success.apply(func=get_from_list, axis=1, args=("MAPQs", "topJ"))
cdr3_success["topV_EditDis"] = cdr3_success.apply(func=get_from_list, axis=1, args=("EditDistance", "topV"))
cdr3_success["topJ_EditDis"] = cdr3_success.apply(func=get_from_list, axis=1, args=("EditDistance", "topJ"))
head_suc = cdr3_success.head()

# %% Plot

def barh(series, title="Title", xlabel="X", ylabel="Y", color='b'):
	plt.figure(figsize=(10, 10))
	counts = series.value_counts(dropna=False)
	plt.barh(counts.index.to_list(), counts, color=color)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#plt.yticks(fontsize=10)
	plt.show()

barh(cdr3_success["topV_Region"], title="SUC topV region", color='g') # Not as distributed, has 5 above 6%.
# TRAV20, TRAV22-1, TRAV25, TRAV

barh(cdr3_error["topV_Region"], title="ERR topV region", color='r')
barh(cdr3_success["topJ_Region"], title="SUC topJ region", color='b') # Pretty distributed
barh(cdr3_error["topJ_Region"], title="ERR topJ region", color='r')

# Count combinations top V and J for error!
# TRAJ28 TRAJ21 ... x TRAV26-1 TRAV26-2 ...


# %%
 
cdr3_positions_file = "/Users/jggatter/Desktop/Projects/SeqWell-TCR/202008192/ref/macFasCDR3bases_v4JG.txt"
fasta = "/Users/jggatter/Desktop/Projects/SeqWell-TCR/20200819/ref/CynoTCRv4.fa"

fasta = pysam.FastaFile(fasta)
cdr3_positions = tcrgo.io.read_cdr3_file(cdr3_positions_file)

for ref in cdr3_positions.keys():
	if 'V' in ref: continue
	if 'C' in ref: continue
	print(ref)
	print(fasta.fetch(ref))
	print(f"{' '*(cdr3_positions[ref]-1)}*")
	print('')

# %%

def hist(series, xlabel="X", ylabel="Y", title="Histogram", bins=20, log=False, density=True, color='g'):
	plt.figure(figsize=(10, 10))
	#counts = series.value_counts(dropna=False)
	plt.hist(series, bins, density=density, log=log, facecolor=color, alpha=0.75)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.grid(True)
	plt.show()

hist(cdr3_success["topV_Score"], title="topV_Score", density=False)
hist(cdr3_success["topJ_Score"], title="topJ_Score", density=False, color='b')
hist(cdr3_success["topV_Phred"], title="topV_Phred", density=False) "A phred score of 20 or above is acceptable"
hist(cdr3_success["topJ_Phred"], title="topJ_Phred", density=False, color='b')
hist(cdr3_success["topV_MAPQ"], title="topV_MAPQ", density=False) "bowtie2 -a Supplementary alignments will also be assigned a MAPQ of 255"
hist(cdr3_success["topJ_MAPQ"], title="topJ_MAPQ", density=False, color='b')
bar(cdr3_success["topV_EditDis"], title="topV_EditDis") # I think filtering out edit distance > 1 or 2 is reasonable
bar(cdr3_success["topJ_EditDis"], title="topJ_EditDis")

# %%
hist(cdr3_error["topV_Score"], title="topV_Score", density=False, color='r')
hist(cdr3_error["topJ_Score"], title="topJ_Score", density=False, color='b')
hist(cdr3_error["topV_Phred"], title="topV_Phred", density=False, color='r')
hist(cdr3_error["topJ_Phred"], title="topJ_Phred", density=False, color='b')
hist(cdr3_error["topV_MAPQ"], title="topV_MAPQ", density=False, color='r')
hist(cdr3_error["topJ_MAPQ"], title="topJ_MAPQ", density=False, color='b')
bar(cdr3_error["topV_EditDis"], title="topV_EditDis")
bar(cdr3_error["topJ_EditDis"], title="topJ_EditDis") # Most failures have AS <80 !

# %% PYSAM

bam = pysam.AlignmentFile("data/190923ShaA_B_raw_bt2_tagged_merged_beadsynth_repaired.bam", 'rb')
bam = pysam.IndexedReads(bam)

for query in queries_cdr3_error.index.to_list():
	read = Read(query)
	read.get_alignments(bam)
	read.print
	
	