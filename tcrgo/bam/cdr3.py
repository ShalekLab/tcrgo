from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from typing import Union, Tuple, List
from collections import defaultdict
import tcrgo.bio as bio

class Transcript(object):
	
	def __init__(self, seq_nt: Seq, query_start_nt: int, \
		V_end_aa: int, J_start_aa: int, orf: int):
		self.seq_nt = seq_nt
		self.seq_aa = None
		self.query_start_nt = query_start_nt
		self.V_end_aa = V_end_aa
		self.J_start_aa = J_start_aa
		self.orf = orf
		self.stops = None

	def get_translation(self):
		end_translation = len(self.seq_nt) - len(self.seq_nt) % 3
		self.seq_aa = self.seq_nt[:end_translation].translate()
	
	def get_stops(self):
		self.stops = bio.find_all(self.seq_aa, '*', 0, len(self.seq_aa))
		return self.stops

	def get_cdr3_starts(self):
		return bio.find_all(
			self.seq_aa, 'C', self.V_end_aa - 8, self.J_start_aa
		)

	def get_cdr3_ends(self):
		return bio.find_all(
			self.seq_aa, 'F', self.V_end_aa, self.J_start_aa + 14
		)

class CDR3(object):
	
	def __init__(self, seq_nt: Seq, start: int, end: int , transcript: Transcript):
		self.seq_nt = seq_nt
		self.seq_aa = None
		self.is_productive = False
		self.start = start
		self.end = end # not end slice
		self.transcript = transcript
		self.score = 0
		
	def get_translation(self):
		if self.end - self.start >= 2:
			end_translation = len(self.seq_nt) - len(self.seq_nt) % 3
			try:
				self.seq_aa = self.seq_nt[:end_translation].translate()
			except TranslationError:
				self.seq_aa = None

	def is_in_frame(self):
		return self.start % 3 == (self.end - 2) % 3

	def score_residues(self):
		if len(self.seq_aa) < 5:
			return -200
		score = 0
		if self.seq_aa[0] == 'C':
			score += 50
		'''
		if self.seq_aa[1] == 'A':
			score += 10
		if self.seq_aa[2] == 'S':
			score += 10
		if self.seq_aa[3] == 'S':
			score += 10
		'''
		if self.seq_aa[-1] == 'F':
			score += 50
		return score

	def calculate_score(self, count_stops: int):
		score = 0
		score += self.score_residues()
		if self.is_in_frame():
			self.is_productive = True
		else:
			score -= 100
		if count_stops > 0:
			self.is_productive = False
		score -= count_stops * 20
		score += len(self.seq_aa)
		self.score = score
