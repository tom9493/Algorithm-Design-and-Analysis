#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.seq1 = ' ' + seq1
		self.seq2 = ' ' + seq2

		if align_length > len(seq1): len1 = len(self.seq1)
		else: len1 = align_length + 1
		if align_length > len(seq2): len2 = len(self.seq2)
		else: len2 = align_length + 1

		# These matrices each take O(nm) space
		self.m = np.full((len2, len1), np.inf)
		self.mb= np.full((len2, len1), "")

		# Banded takes O(kn) time
		if banded:
			if abs(len1-len2) > 3:
				score, alignment1, alignment2 = np.inf, "No Alignment Available", "No Alignment Available"
			else: score, alignment1, alignment2 = self.bandedAlign(len1, len2)

		# unBanded takes O(nm) time
		else:
			score, alignment1, alignment2 = self.nonBandedAlign(len1, len2)

		alignment1 = alignment1[0:100]
		alignment2 = alignment2[0:100]
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	# Aligns the gene sequences with a certain bandwidth
	# O(kn) time and O(1) space as the matrices were already made
	def bandedAlign(self, len1, len2):
		d = 3
		A = min(len1, len2)

		# First loop is length of shorter sequence. Second is k = 2d + 1
		# Constant work for each loop so O(kn)
		for i in range(0, A):
			for j in range(0, d+1):	# Does 3 cells to the right of diagonal
				if i+j < len1:
					if self.seq2[i] == self.seq1[i+j]: match = -3
					else: match = 1

					a = self.m[i][i+j-1] + 5
					b = self.m[i-1][i+j-1] + match
					c = self.m[i-1][i+j] + 5

					e = min(a, b, c)
					self.m[i][i+j] = e

					if e == a: self.mb[i][i+j] = "a"
					elif e == c: self.mb[i][i+j] = "c"
					else: self.mb[i][i+j] = "b"

					if i == 0 and j == 0: self.m[i][i] = 0
					if i == 0 and j != 0:
						self.m[i][j] = self.m[0][j-1] + 5
						self.mb[i][j] = "a"

				else: continue

			for j in range(1, d+1): # Does 3 cells below the diagonal
				if i+j < len2:
					if self.seq2[i+j] == self.seq1[i]: match = -3
					else: match = 1

					a = self.m[i+j][i-1] + 5
					b = self.m[i+j-1][i-1] + match
					c = self.m[i+j-1][i] + 5

					e = min(a, b, c)
					self.m[i+j][i] = e

					if e == a: self.mb[i+j][i] = "a"
					elif e == c: self.mb[i+j][i] = "c"
					else: self.mb[i+j][i] = "b"

					if i != 0 and j == 0 and i < d+1:
						self.m[i][j] = self.m[i-1][0] + 5
						self.mb[i][j] = "c"

				else: continue

		score = self.m[len2 - 1][len1 - 1]
		alignment1, alignment2 = self.getAlignments(len1, len2)
		return score, alignment1, alignment2

	# Aligns the gene sequences by filling out entire alignment table
	# O(nm) time and O(1) space
	def nonBandedAlign(self, len1, len2):
		for i in range(len2):
			for j in range(len1):
				if self.seq1[j] == self.seq2[i]: match = -3
				else: match = 1

				a = self.m[i][j-1] + 5
				b = self.m[i-1][j-1] + match
				c = self.m[i-1][j] + 5

				e = min(a, b, c)
				self.m[i][j] = e

				if e == a: self.mb[i][j] = "a"
				elif e == c: self.mb[i][j] = "c"
				else: self.mb[i][j] = "b"

				if i == 0 and j == 0: self.m[i][j] = 0
				if i != 0 and j == 0:
					self.m[i][j] = self.m[i-1][0] + 5
					self.mb[i][j] = "c"
				if i == 0 and j != 0:
					self.m[i][j] = self.m[0][j-1] + 5
					self.mb[i][j] = "a"

		score = self.m[len2-1][len1-1]
		alignment1, alignment2 = self.getAlignments(len1, len2)
		return score, alignment1, alignment2

	# Traces back through the matrix to find the alignment
	# Worst case O(n+m) time, O(2n+1) for banded, O(1) space
	def getAlignments(self, len1, len2):
		i, j = len2-1, len1-1
		a1 = ""
		a2 = ""

		while i != 0 and j != 0:
			if self.mb[i][j] == "a":
				a1 = self.seq1[j] + a1
				a2 = "-" + a2
				j -= 1
			if self.mb[i][j] == "c":
				a1 = "-" + a1
				a2 = self.seq2[i] + a2
				i -= 1
			else:
				a1 = self.seq1[j] + a1
				a2 = self.seq2[i] + a2
				j -= 1
				i -= 1

		return a1, a2