#! coding: utf-8

import os

from argparse import ArgumentParser
from Bio.Align import PairwiseAligner
from collections import defaultdict
from collections import namedtuple

__author__ = "wei"
__date__   = "20210521"
__email__  = "hanwei@shanghaitech.edu.cn"


class Comdline:
	"""
	parse comdline argparse
	"""
	def __init__(self):
		self.args = {}
	
	def argparse(self):
		comdline = ArgumentParser()
		comdline.add_argument('-i', '--input', type=str,
			help='input filename in fasta format, required')
		comdline.add_argument('-o', '--output', type=str,
			help='output filename, required')
		comdline.add_argument('-c', '--cutoff', type=float,
			default=0.9,
			help='sequence identity threshold, default 0.9')
		self.args = comdline.parse_args()


class CDHIT:
	"""
	According to the CDHIT paper, try using Python to implatement the tool.
	This is just to help me to understand it, and to salute CDHIT.
	"""
	def __init__(self, seqfile, thre, blosum, prefix, kmer=5):
		self.seqfile = seqfile
		self.thre = thre # cluster threshold
		self.seqlist = []
		self.prefix = prefix
		self.blosum = blosum
		self.kmer = kmer


	def readfasta(self):
		"""
		open and read sequence file
		Result:
			type -> list
			data -> [(SequenceName, SequeneLength, Sequence), ...]
		"""
		with open(self.seqfile) as seqf:
			seqs = seqf.read().split('>')[1:]
		for seq in seqs:
			lines = seq.split('\n')
			seqname = lines[0]
			seqaads = ''.join(lines[1:])
			seqleng = len(seqaads)
			self.seqlist.append((seqname, seqleng, seqaads))
	
	
	def sort_seq_length(self):
		"""
		sort from smallest to largest
		"""
		self.seqlist = sorted(
			self.seqlist,
			key=lambda x: x[1],
			reverse=False
		)


	def kmercount(self, seq):
		"""
		word count
		split sequence to word by kmer, and count
		"""
		wordcount = defaultdict(int)
		seqlen = len(seq)
		for i in range(seqlen-self.kmer):
			word = seq[i:i+self.kmer] # split sequence
			wordcount[word] += 1
		return wordcount


	def short_filter(self, refseq, seq, reflen, th_val):
		"""
		short word filter
		qulickly judhe whether the clustering requirements are met
		"""
		refcount = self.kmercount(refseq)
		seqcount = self.kmercount(seq)
		overcount = 0
		# determine sequence length to reduce time complexity
		if len(refcount) < len(seqcount):
			for mer, c1 in refcount.items():
				c2 = seqcount.get(mer, 0)
				overcount += min(c1, c2)
		else:
			for mer, c1 in seqcount.items():
				c2 = refcount.get(mer, 0)
				overcount += min(c1, c2)
		return False if overcount < th_val else True


	def calculate_ident(self, refseq, seq):
		#align = PairwiseAlign(refseq, seq, self.blosum)
		#align.align()
		#res = align.recover_alignments()
		aligner = PairwiseAligner()
		alignments = aligner.align(refseq, seq)
		aligned = alignments[0].aligned
		query = alignments[0].query
		target = alignments[0].target
		ident = sum(j-i for i, j in aligned[0])
		ident = ident / min(len(query), len(target))
		return ident


	def cluster(self):
		"""
		clustering according to dentity of the sequence
		steps:
			1: sort by length of sequence  *
			2: short word filtering        ***
			3: identity filter             **
		"""
		allcluster, cluseq = [], []
		while self.seqlist:
			refname, reflen, refseq = self.seqlist.pop()
			cluseq.append('>'+refname+'\n'+refseq+'\n')
			
			# transform from identity to precentage difference
			# if identity = 0.9, then difference=1-0.9
			diff = 1 - self.thre
			# short word filtering threshold
			th_val = int(reflen * (1 -self.kmer*diff))

			overlength = len(self.seqlist)
			clusterindex = []
			cluster = [(refname, reflen, "*", "refer")]
			for i in range(overlength-1, -1, -1):
				name, leng, seq = self.seqlist[i]
				# shor word filter
				ismatch = self.short_filter(refseq, seq, reflen, th_val)
				if not ismatch: continue
				# calculate identity
				ident = self.calculate_ident(refseq, seq)
				if ident < self.thre: continue
				clusterindex.append(i)
				cluster.append((name, leng, ident, 'cluster'))
			for i in clusterindex:
				del self.seqlist[i]
			allcluster.append(cluster)
		return allcluster, cluseq

	def writer_cluster(self, cluster):
		"""
		write cluster information to file
		"""
		cluout = open(self.prefix+'.cluster', 'w')
		count = 0
		for clu in cluster:
			head = '>Cluster ' + str(count) + '\n'
			cluout.write(head)
			count += 1
			num = 0
			for name, leng, ident, _type in clu:
				if _type == "cluster":
					line = "{: <4d}{:<5d} >{:<50s}... at {:.2f}%\n".\
						format(num, leng, name, ident*100)
				else:
					line = "{: <4d}{:<5d} >{:<50s}... *\n".\
						format(num, leng, name)
				num += 1
				cluout.write(line)
		cluout.close()

	def writer_cluseq(self, cluseq):
		"""
		write refersence sequence to file
		"""
		seqout = open(self.prefix+'.fasta', 'w')
		seqout.writelines(cluseq)
		seqout.close()
		


class PairwiseAlign:
	"""
	pairwise sequence alignment
	"""
	def __init__(
		self,
		seqA,
		seqB,
		blosum,
		_open=-2,
		_extend=-1
	):
		self.seqA = seqA
		self.seqB = seqB
		self.blosum = blosum
		self.open = _open
		self.extend = _extend
		self.ident = 0
		self.alignres = {}

	def align(self):
		lenA, lenB = len(self.seqA), len(self.seqB)
		score_matrix, trace_matrix = [], []
		# Initialize score_matrix and trace_matrix
		for _ in range(lenA+1):
			score_matrix.append([None] * (lenB + 1))
			trace_matrix.append([None] * (lenB + 1))
		# Initialize first row and column with gap scores
		for i in range(lenA + 1):
			score = self.gap_score(i)
			score_matrix[i][0] = score
		for i in range(lenB + 1):
			score = self.gap_score(i)
			score_matrix[0][i] = score
		# Fill in the score matrix
		for row in range(1, lenA + 1):
			for col in range(1, lenB + 1):
				nogap_score = (
					score_matrix[row-1][col-1]
					+ self.blosum[self.seqA[row-1] + self.seqB[col-1]]
				)
				row_open = score_matrix[row][col-1] + self.open
				row_extend = max(
					score_matrix[row][x] + self.gap_score(col-x) for x in range(col)
				)
				col_open = score_matrix[row-1][col] + self.open
				col_extend = max(
					score_matrix[x][col] + self.gap_score(row-x) for x in range(row)
				)
				best = max(nogap_score, row_open, row_extend, col_open, col_extend)
				score_matrix[row][col] = best
				trace_score = 0
				if self.rint(nogap_score) == self.rint(best):
					trace_score += 2
				if self.rint(row_open) == self.rint(best):
					trace_score += 1
				if self.rint(row_extend) == self.rint(best):
					trace_score += 8
				if self.rint(col_open) == self.rint(best):
					trace_score += 4
				if self.rint(col_extend) == self.rint(best):
					trace_score += 16
				trace_matrix[row][col] = trace_score

		self.alignres = {
			'score_matrix': score_matrix,
			'trace_matrix': trace_matrix,
			'best_score': best,
		}

	
	def recover_alignments(self, reverse=False):
		lenA, lenB = len(self.seqA), len(self.seqB)
		ali_seqA, ali_seqB = self.seqA[0:0], self.seqB[0:0]
		gap_char = '-'
		tracebacks = []
		in_process = []
		
		score_matrix = self.alignres['score_matrix']
		trace_matrix = self.alignres['trace_matrix']
		score = self.alignres['best_score']
		nrows, ncols = len(score_matrix), len(score_matrix[0])
		row = nrows - 1
		col = ncols - 1
		begin = 0
		end = None
		in_process += [
			(ali_seqA, ali_seqB, end, row, col, False, trace_matrix[row][col])
		]
		while in_process:
			dead_end = False
			ali_seqA, ali_seqB, end, row, col, col_gap, trace = in_process.pop()
	
			while (row > 0 or col > 0) and not dead_end:
				cache = (ali_seqA[:], ali_seqB[:], end, row, col, col_gap)
	
				# If trace is empty we have reached at least one border of the
				# matrix or the end of a local alignment. Just add the rest of
				# the sequence(s) and fill with gaps if necessary.
				if not trace:
					if col and col_gap:
						dead_end = True
					else:
						ali_seqA, ali_seqB = self._finish_backtrace(
							self.seqA, self.seqB, ali_seqA, ali_seqB, row, col, gap_char
						)
					break
				elif trace % 2 == 1:  # = row open = open gap in seqA
					trace -= 1
					if col_gap:
						dead_end = True
					else:
						col -= 1
						ali_seqA += gap_char
						ali_seqB += self.seqB[col : col + 1]
						col_gap = False
				elif trace % 4 == 2:  # = match/mismatch of seqA with seqB
					trace -= 2
					row -= 1
					col -= 1
					ali_seqA += self.seqA[row : row + 1]
					ali_seqB += self.seqB[col : col + 1]
					col_gap = False
				elif trace % 8 == 4:  # = col open = open gap in seqB
					trace -= 4
					row -= 1
					ali_seqA += self.seqA[row : row + 1]
					ali_seqB += gap_char
					col_gap = True
				elif trace in (8, 24):  # = row extend = extend gap in seqA
					trace -= 8
					if col_gap:
						dead_end = True
					else:
						col_gap = False
						# We need to find the starting point of the extended gap
						x = self._find_gap_open(
							self.seqA,
							self.seqB,
							ali_seqA,
							ali_seqB,
							end,
							row,
							col,
							col_gap,
							gap_char,
							score_matrix,
							trace_matrix,
							in_process,
							col,
							row,
							"col",
							score,
						)
						ali_seqA, ali_seqB, row, col, in_process, dead_end = x
				elif trace == 16:  # = col extend = extend gap in seqB
					trace -= 16
					col_gap = True
					x = self._find_gap_open(
						self.seqA,
						self.seqB,
						ali_seqA,
						ali_seqB,
						end,
						row,
						col,
						col_gap,
						gap_char,
						score_matrix,
						trace_matrix,
						in_process,
						row,
						col,
						"row",
						score,
					)
					ali_seqA, ali_seqB, row, col, in_process, dead_end = x
	
				if trace:  # There is another path to follow...
					cache += (trace,)
					in_process.append(cache)
				trace = trace_matrix[row][col]
			if not dead_end:
				if not reverse:
					tracebacks.append((ali_seqA[::-1], ali_seqB[::-1], score, begin, end))
				else:
					tracebacks.append((ali_seqB[::-1], ali_seqA[::-1], score, begin, end))
		return self._clean_alignments(tracebacks)

	def gap_score(self, gaplen):
		score = self.open + gaplen * self.extend
		return score


	def rint(self, x, precision=1000):
		return int(x * precision + 0.5)


	def _finish_backtrace(
		self,
		sequenceA,
		sequenceB,
		ali_seqA,
		ali_seqB,
		row,
		col,
		gap_char
	):
	    """Add remaining sequences and fill with gaps if necessary (PRIVATE)."""
	    if row:
	        ali_seqA += sequenceA[row - 1 :: -1]
	    if col:
	        ali_seqB += sequenceB[col - 1 :: -1]
	    if row > col:
	        ali_seqB += gap_char * (len(ali_seqA) - len(ali_seqB))
	    elif col > row:
	        ali_seqA += gap_char * (len(ali_seqB) - len(ali_seqA))
	    return ali_seqA, ali_seqB


	def _clean_alignments(self, alignments):
	    """Take a list of alignments and return a cleaned version (PRIVATE).
	
	    Remove duplicates, make sure begin and end are set correctly, remove
	    empty alignments.
	    """
	    Alignment = namedtuple("Alignment", ("seqA, seqB, score, start, end"))
	    unique_alignments = []
	    for align in alignments:
	        if align not in unique_alignments:
	            unique_alignments.append(align)
	    i = 0
	    while i < len(unique_alignments):
	        seqA, seqB, score, begin, end = unique_alignments[i]
	        # Make sure end is set reasonably.
	        if end is None:  # global alignment
	            end = len(seqA)
	        elif end < 0:
	            end = end + len(seqA)
	        # If there's no alignment here, get rid of it.
	        if begin >= end:
	            del unique_alignments[i]
	            continue
	        unique_alignments[i] = Alignment(seqA, seqB, score, begin, end)
	        i += 1
	    return unique_alignments


	def _find_gap_open(
		self,
		sequenceA,
		sequenceB,
		ali_seqA,
		ali_seqB,
		end,
		row,
		col,
		col_gap,
		gap_char,
		score_matrix,
		trace_matrix,
		in_process,
		target,
		index,
		direction,
		best_score,
	):
		"""Find the starting point(s) of the extended gap (PRIVATE)."""
		dead_end = False
		target_score = score_matrix[row][col]
		for n in range(target):
			if direction == "col":
				col -= 1
				ali_seqA += gap_char
				ali_seqB += sequenceB[col : col + 1]
			else:
				row -= 1
				ali_seqA += sequenceA[row : row + 1]
				ali_seqB += gap_char
			actual_score = score_matrix[row][col] + self.gap_score(n + 1)
			if self.rint(actual_score) == self.rint(target_score) and n > 0:
				if not trace_matrix[row][col]:
					break
				else:
					in_process.append(
						(
							ali_seqA[:],
							ali_seqB[:],
							end,
							row,
							col,
							col_gap,
							trace_matrix[row][col],
						)
					)
			if not trace_matrix[row][col]:
				dead_end = True
		return ali_seqA, ali_seqB, row, col, in_process, dead_end
	


BLOSUM62 = {
	'AA': 4, 'AC': 0, 'AD': -2, 'AE': -1, 'AF': -2, 'AG': 0, 
	'AH': -2, 'AI': -1, 'AK': -1, 'AL': -1, 'AM': -1, 'AN': -2, 
	'AP': -1, 'AQ': -1, 'AR': -1, 'AS': 1, 'AT': 0, 'AV': 0, 
	'AW': -3, 'AX': 0, 'AY': -2, 'CA': 0, 'CC': 9, 'CD': -3, 
	'CE': -4, 'CF': -2, 'CG': -3, 'CH': -3, 'CI': -1, 'CK': -3, 
	'CL': -1, 'CM': -1, 'CN': -3, 'CP': -3, 'CQ': -3, 'CR': -3, 
	'CS': -1, 'CT': -1, 'CV': -1, 'CW': -2, 'CX': -2, 'CY': -2, 
	'DA': -2, 'DC': -3, 'DD': 6, 'DE': 2, 'DF': -3, 'DG': -1, 
	'DH': -1, 'DI': -3, 'DK': -1, 'DL': -4, 'DM': -3, 'DN': 1, 
	'DP': -1, 'DQ': 0, 'DR': -2, 'DS': 0, 'DT': -1, 'DV': -3, 
	'DW': -4, 'DX': -1, 'DY': -3, 'EA': -1, 'EC': -4, 'ED': 2, 
	'EE': 5, 'EF': -3, 'EG': -2, 'EH': 0, 'EI': -3, 'EK': 1, 
	'EL': -3, 'EM': -2, 'EN': 0, 'EP': -1, 'EQ': 2, 'ER': 0, 
	'ES': 0, 'ET': -1, 'EV': -2, 'EW': -3, 'EX': -1, 'EY': -2, 
	'FA': -2, 'FC': -2, 'FD': -3, 'FE': -3, 'FF': 6, 'FG': -3, 
	'FH': -1, 'FI': 0, 'FK': -3, 'FL': 0, 'FM': 0, 'FN': -3, 
	'FP': -4, 'FQ': -3, 'FR': -3, 'FS': -2, 'FT': -2, 'FV': -1, 
	'FW': 1, 'FX': -1, 'FY': 3, 'GA': 0, 'GC': -3, 'GD': -1, 
	'GE': -2, 'GF': -3, 'GG': 6, 'GH': -2, 'GI': -4, 'GK': -2, 
	'GL': -4, 'GM': -3, 'GN': 0, 'GP': -2, 'GQ': -2, 'GR': -2, 
	'GS': 0, 'GT': -2, 'GV': -3, 'GW': -2, 'GX': -1, 'GY': -3, 
	'HA': -2, 'HC': -3, 'HD': -1, 'HE': 0, 'HF': -1, 'HG': -2, 
	'HH': 8, 'HI': -3, 'HK': -1, 'HL': -3, 'HM': -2, 'HN': 1, 
	'HP': -2, 'HQ': 0, 'HR': 0, 'HS': -1, 'HT': -2, 'HV': -3, 
	'HW': -2, 'HX': -1, 'HY': 2, 'IA': -1, 'IC': -1, 'ID': -3, 
	'IE': -3, 'IF': 0, 'IG': -4, 'IH': -3, 'II': 4, 'IK': -3, 
	'IL': 2, 'IM': 1, 'IN': -3, 'IP': -3, 'IQ': -3, 'IR': -3, 
	'IS': -2, 'IT': -1, 'IV': 3, 'IW': -3, 'IX': -1, 'IY': -1, 
	'KA': -1, 'KC': -3, 'KD': -1, 'KE': 1, 'KF': -3, 'KG': -2, 
	'KH': -1, 'KI': -3, 'KK': 5, 'KL': -2, 'KM': -1, 'KN': 0, 
	'KP': -1, 'KQ': 1, 'KR': 2, 'KS': 0, 'KT': -1, 'KV': -2, 
	'KW': -3, 'KX': -1, 'KY': -2, 'LA': -1, 'LC': -1, 'LD': -4, 
	'LE': -3, 'LF': 0, 'LG': -4, 'LH': -3, 'LI': 2, 'LK': -2, 
	'LL': 4, 'LM': 2, 'LN': -3, 'LP': -3, 'LQ': -2, 'LR': -2, 
	'LS': -2, 'LT': -1, 'LV': 1, 'LW': -2, 'LX': -1, 'LY': -1, 
	'MA': -1, 'MC': -1, 'MD': -3, 'ME': -2, 'MF': 0, 'MG': -3, 
	'MH': -2, 'MI': 1, 'MK': -1, 'ML': 2, 'MM': 5, 'MN': -2, 
	'MP': -2, 'MQ': 0, 'MR': -1, 'MS': -1, 'MT': -1, 'MV': 1, 
	'MW': -1, 'MX': -1, 'MY': -1, 'NA': -2, 'NC': -3, 'ND': 1, 
	'NE': 0, 'NF': -3, 'NG': 0, 'NH': 1, 'NI': -3, 'NK': 0, 
	'NL': -3, 'NM': -2, 'NN': 6, 'NP': -2, 'NQ': 0, 'NR': 0, 
	'NS': 1, 'NT': 0, 'NV': -3, 'NW': -4, 'NX': -1, 'NY': -2, 
	'PA': -1, 'PC': -3, 'PD': -1, 'PE': -1, 'PF': -4, 'PG': -2, 
	'PH': -2, 'PI': -3, 'PK': -1, 'PL': -3, 'PM': -2, 'PN': -2, 
	'PP': 7, 'PQ': -1, 'PR': -2, 'PS': -1, 'PT': -1, 'PV': -2, 
	'PW': -4, 'PX': -2, 'PY': -3, 'QA': -1, 'QC': -3, 'QD': 0, 
	'QE': 2, 'QF': -3, 'QG': -2, 'QH': 0, 'QI': -3, 'QK': 1, 
	'QL': -2, 'QM': 0, 'QN': 0, 'QP': -1, 'QQ': 5, 'QR': 1, 
	'QS': 0, 'QT': -1, 'QV': -2, 'QW': -2, 'QX': -1, 'QY': -1, 
	'RA': -1, 'RC': -3, 'RD': -2, 'RE': 0, 'RF': -3, 'RG': -2, 
	'RH': 0, 'RI': -3, 'RK': 2, 'RL': -2, 'RM': -1, 'RN': 0, 
	'RP': -2, 'RQ': 1, 'RR': 5, 'RS': -1, 'RT': -1, 'RV': -3, 
	'RW': -3, 'RX': -1, 'RY': -2, 'SA': 1, 'SC': -1, 'SD': 0, 
	'SE': 0, 'SF': -2, 'SG': 0, 'SH': -1, 'SI': -2, 'SK': 0, 
	'SL': -2, 'SM': -1, 'SN': 1, 'SP': -1, 'SQ': 0, 'SR': -1, 
	'SS': 4, 'ST': 1, 'SV': -2, 'SW': -3, 'SX': 0, 'SY': -2, 
	'TA': 0, 'TC': -1, 'TD': -1, 'TE': -1, 'TF': -2, 'TG': -2, 
	'TH': -2, 'TI': -1, 'TK': -1, 'TL': -1, 'TM': -1, 'TN': 0, 
	'TP': -1, 'TQ': -1, 'TR': -1, 'TS': 1, 'TT': 5, 'TV': 0, 
	'TW': -2, 'TX': 0, 'TY': -2, 'VA': 0, 'VC': -1, 'VD': -3, 
	'VE': -2, 'VF': -1, 'VG': -3, 'VH': -3, 'VI': 3, 'VK': -2, 
	'VL': 1, 'VM': 1, 'VN': -3, 'VP': -2, 'VQ': -2, 'VR': -3, 
	'VS': -2, 'VT': 0, 'VV': 4, 'VW': -3, 'VX': -1, 'VY': -1, 
	'WA': -3, 'WC': -2, 'WD': -4, 'WE': -3, 'WF': 1, 'WG': -2, 
	'WH': -2, 'WI': -3, 'WK': -3, 'WL': -2, 'WM': -1, 'WN': -4, 
	'WP': -4, 'WQ': -2, 'WR': -3, 'WS': -3, 'WT': -2, 'WV': -3, 
	'WW': 11, 'WX': -2, 'WY': 2, 'XA': 0, 'XC': -2, 'XD': -1, 
	'XE': -1, 'XF': -1, 'XG': -1, 'XH': -1, 'XI': -1, 'XK': -1, 
	'XL': -1, 'XM': -1, 'XN': -1, 'XP': -2, 'XQ': -1, 'XR': -1, 
	'XS': 0, 'XT': 0, 'XV': -1, 'XW': -2, 'XX': -1, 'XY': -1, 
	'YA': -2, 'YC': -2, 'YD': -3, 'YE': -2, 'YF': 3, 'YG': -3, 
	'YH': 2, 'YI': -1, 'YK': -2, 'YL': -1, 'YM': -1, 'YN': -2, 
	'YP': -3, 'YQ': -1, 'YR': -2, 'YS': -2, 'YT': -2, 'YV': -1, 
	'YW': 2, 'YX': -1, 'YY': 7,
}


def main():
	# parse comdline
	comdline = Comdline()
	comdline.argparse()
	args = comdline.args
	infile = args.input
	outfile = args.output
	cutoff = args.cutoff
	# check comdline
	if (not infile) or (not outfile):
		os.system('python cdhit.py -h')
		print('\033[31;1minput or output arguments are missing.\033[0m')
		os._exit(0)
	# run cd-hit	
	cdhit = CDHIT(infile, cutoff, BLOSUM62, outfile)
	cdhit.readfasta()
	cdhit.sort_seq_length()
	cluster, cluseq = cdhit.cluster()
	cdhit.writer_cluster(cluster)
	cdhit.writer_cluseq(cluseq)
	
if __name__ == "__main__":
	main()
