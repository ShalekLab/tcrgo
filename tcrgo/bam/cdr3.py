"""
An idea for later possibly. A lot of methods involving CDR3s also alignment data so it seems 
just about all methods would remain in read.py. Advantage of this is that less memory would 
be allocated as non-CDR3 reads wouldn't take up these fields
"""

class CDR3(object):
	
	def __init__(self):
		self.sequence_nucleotide = None
		self.sequence_aminoacid = None
		self.is_complete = False
		self.ref_start = None
		self.query_start = None
		self.ref_end = None
		self.query_end = None
		self.ref_seq_J = "" #?
		self.ref_seq_V = "" #?