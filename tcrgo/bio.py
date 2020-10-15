def levenshtein_distance(sequence1, sequence2):
	edit_distances = []
	# Create matrix
	for i in range(len(sequence1)+1):
		edit_distances.append([0] * len(sequence2)+1)
	# Initialize ED(a, empty string prefix of b) = i
	for i in range(len(sequence1)+1):
		edit_distances[i][0] = i
	# Initialize ED(b, empty string prefix of a) = j
	for j in range(len(sequence2)+1):
		edit_distances[0][j] = j

	for i in range(1, len(sequence1)+1):
		for j in range(1, len(sequence2)+1):
			distance_horizontal = edit_distances[i][j-1] + 1 # indel from b to a, penalty +1
			distance_vertical = edit_distances[i-1][j] + 1 # indel from a to b, penalty +1
			# substitution at end of a and b
			if sequence1[i-1] == sequence2[j-1]: # If characters match, no penalty
				distance_diagonal = edit_distances[i-1][j-1]
			else: # Characters do not match, penalty +1
				distance_diagonal = edit_distances[i-1][j-1] + 1
			edit_distances[i][j] = min(distance_horizontal, distance_vertical, distance_diagonal)
	return edit_distances[-1][-1] # Bottom right corner is the overall edit distance