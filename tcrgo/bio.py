from Bio.Seq import Seq
from typing import List, Optional, Tuple, Dict, Set, Any, Union
import pysam
from collections import Counter, defaultdict
import numpy as np
Array = np.ndarray
AlignedSegment = pysam.libcalignedsegment.AlignedSegment

# Faster solution. str is faster than Seq
# https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def levenshtein_distance(seq1: Union[str, Seq], seq2: Union[str, Seq]):
	# Potentially flip the problem to ensure the longer sequence is seq1
	if len(seq1) < len(seq2):
		return levenshtein_distance(seq2, seq1)
	# If one of the sequences is 0, return the length of the other (potentially also zero)
	if len(seq2) == 0:
		return len(seq1)

	previous_row = range(len(seq2) + 1)
	for i, c1 in enumerate(seq1):
		current_row = [i + 1]
		for j, c2 in enumerate(seq2):
			# j+1 instead of j since previous_row and current_row are one character longer than seq2
			insertions = previous_row[j + 1] + 1
			deletions = current_row[j] + 1
			substitutions = previous_row[j] + (c1 != c2) # If characters match, no penalty. Branchless?
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row

	return previous_row[-1]

# Slightly slower. From Jon Hopkins online course.
# I found str types to be 30-40% faster than Seq.
def _OLD_levenshtein_distance(sequence1: Union[str, Seq], sequence2: Union[str, Seq]) -> int:
	edit_distances = []
	# Create matrix
	for i in range(len(sequence1)+1):
		edit_distances.append([0] * (len(sequence2)+1))
	# Initialize ED(a, empty string prefix of b) = i
	for i in range(len(sequence1)+1):
		edit_distances[i][0] = i
	# Initialize ED(b, empty string prefix of a) = j
	for j in range(len(sequence2)+1):
		edit_distances[0][j] = j
	for i, c1 in enumerate(sequence1, start=1):
		for j, c2 in enumerate(sequence2, start=1):
			distance_horizontal = edit_distances[i][j-1] + 1 # indel from b to a, penalty +1
			distance_vertical = edit_distances[i-1][j] + 1 # indel from a to b, penalty +1
			distance_diagonal = edit_distances[i-1][j-1] + (c1 != c2)
			edit_distances[i][j] = min(distance_horizontal, distance_vertical, distance_diagonal)
	return edit_distances[-1][-1] # Bottom right corner is the overall edit distance

def hamming_distance(sequence1: Union[str, Seq], sequence2: Union[str, Seq]) -> int:
	if len(sequence1) != len(sequence2):
		raise ValueError("Sequences are not of the same length.")
	count_mismatches = 0
	for i in range(len(sequence1)):
		if sequence1[i] != sequence2[i]:
			count_mismatches += 1
	return count_mismatches

def levenshtein_distance_bcumi(seq1:str, seq2:str) -> bool:
	# Check if one but not other ends with T. If so, check for 1 deletion.
	# Else check for 1 ins?
	# else <= 1 sub.
	# Find the best order of doing these things.
	pass


def get_frame_translations(sequence_nt: Seq) -> Tuple[Seq, ...]:
	sequences_aa = [None] * 3
	for frame in range(3):
		end = len(sequence_nt) - ((len(sequence_nt) - frame) % 3)
		sequences_aa[frame] = sequence_nt[frame:end].translate()
	return tuple(sequences_aa)

def find_all(seq: Seq, residue: str, \
	start_index: int, end_slice: int) -> Tuple[Optional[int], ...]:
	"""Search for a residue and return a list of positions for each frame"""
	positions = list()
	start = start_index
	end = end_slice
	while start < end:
		position = seq.find(residue, start, end)
		if position == -1:
			break
		positions.append(position)
		start = position + 1
	return tuple(positions)

def find_all_in_frames(seqs_aa: Tuple[Seq, ...], residue: str, \
	start_index: int=0, end_slice: int=None, as_aminoacid_positions: bool=False, \
	are_from_end: bool=False) -> Tuple[Tuple[Optional[int], ...], ...]:
	""" """
	mult = 3 >> as_aminoacid_positions # 1 if True, 3 if False
	has_frame_offset = not as_aminoacid_positions
	frames_positions = [None] * 3
	for f in range(3):
		seq = seqs_aa[f]
		if end_slice is None:
			end_slice = len(seq)
		elif are_from_end:
			end_slice = len(seq) - end_slice
		if are_from_end:
			start_index = len(seq) - start_index
		positions = find_all(seq, residue, start_index, end_slice)
		frames_positions[f] = tuple(
			[pos * mult + has_frame_offset * f for pos in positions]
		)
	return tuple(frames_positions)

def get_frame_max(frames: Tuple[Tuple[Any, ...], ...]) -> Set[int]:
	frames_max = set()
	count_max=0
	for i in range(len(frames)):
		count = len(frames[i])
		if count > count_max:
			count_max = count
			frames_max = set([i])
		elif count == count_max:
			frames_max.add(i)
	return frames_max

def get_frame_min(frames: Tuple[Tuple[Any, ...], ...]) -> Set[int]:
	frames_min = set([0])
	count_min=len(frames[0])
	for i in range(1, len(frames)):
		count = len(frames[i])
		if count < count_min:
			count_min = count
			frames_min = set([i])
		elif count == count_min:
			frames_min.add(i)
	return frames_min

def contains_indels(alignment: AlignedSegment) -> bool:
	"""Return true if insertions or deletions present in CIGAR"""
	return alignment.get_cigar_stats()[0][1] > 0 \
		or alignment.get_cigar_stats()[0][2] > 0

def resolve_ambiguous_base(bases: List[str]) -> str:
	""""
	Based on the tied bases, resolve according to the index table
	N and '_' (No base) automatically lose out.
	"""
	if 'N' in bases:
		bases.remove('N')
	if '_' in bases:
		bases.remove('_')
	if len(bases) == 1:
		return bases[0]
	ambiguous = (
		'X', 'A', 'C', 'M', 
		'G', 'R', 'S', 'V', 
		'T', 'W', 'Y', 'H',
		'K', 'D', 'B', 'X'
	)
	index = 0
	nucleotides = ('A', 'C', 'G', 'T')
	for i in range(len(nucleotides)):
		if nucleotides[i] in bases:
			index += 2**i
	return ambiguous[index]	

def get_consensus_sequence(seq_counts: Dict[Seq, int]) -> Seq:
	"""
	Assuming the sequences are aligned at the beginning,
	Produce a consensus sequence. Greedily lengthens and provides
	ambiguous characters if necessary.
	"""
	consensus = []
	lengths = [len(seq) for seq in seq_counts.keys()]
	for i in range(max(lengths)):
		bases = Counter()
		for seq, count in seq_counts.items():
			if i < len(seq):
				bases[seq[i]] += count
			else:
				bases['_'] += count
		top_count = max(bases.values())
		bases = [base for base,count in bases.items() if count == top_count]
		if len(bases) == 1:
			base = bases[0]
		else:
			base = resolve_ambiguous_base(bases)
		if base != '_':
			consensus.append(base)
		else:
			break
	return Seq(''.join(consensus))

def incorporate_mismatches_deletions(md_string: str, output: List[str], start: int=0) -> List[str]:
	position = start
	num = ""
	for char in md_string:
		if  47 < ord(char) < 58: # If number 0-9
			num += char
			is_deletion = False
		else:
			if num != '':
				end = position + int(num)
				for i in range(position, end):
					output[i] = '|'
				position = end
			if char == '^':
				is_deletion = True
			elif is_deletion:
				output[position] = 'D'
				position += 1
			else:
				output[position] = ' '
				position += 1
			num = ""
	if num != '':
		for i in range(position, position + int(num)):
			output[i] = '|'
	return output

# TODO: Handle insertions and deletions properly.
# TODO: displayed reference sequence will need to be adjusted visually
# in another method probably.
# Insertions add space in reference
# deletions remove bases. these bases could either be displayed in cigar
# or above the reference sequence.
def visualize_alignment(alignment: AlignedSegment) -> str:
	output = ['?'] * (alignment.infer_read_length() \
		+ alignment.get_cigar_stats()[0][1] \
		+ alignment.get_cigar_stats()[0][2])
	num = ""
	position = 0
	read_M = False
	for char in alignment.cigarstring:
		if  47 < ord(char) < 58: # If number 0-9
			num += char
		else:
			if char == 'S':
				for i in range(position, position + int(num)):
					output[i] = '-'
			elif char == 'M' and not read_M:
				incorporate_mismatches_deletions(alignment.get_tag("MD"), output, position)
				read_M = True
			elif char == 'I':
				for i in range(position, position + int(num)):
					output[i] = 'I'
			position += int(num)
			num = ""	
	return ''.join(output)

def visualize_aminoacids(sequence: Seq) -> str:
	output = [' '] * (len(sequence) * 3)
	for i in range(len(sequence)):
		output[i*3] = sequence[i]
	return ''.join(output)
