#a dynamic programming algorithm that finds the optimal global alignment of two input DNA sequences. 
# The program will be run from the command-line as:
# <programfilename> <fastafilename1> <fastafilename2> <substitutionmatrixfile> 
# <gap-penalty>
# ex: python task1.py S2.txt S3.txt subs.txt -500

import math
import numpy as np
from sklearn import decomposition
import wave, struct
import pickle
from pprint import pprint
import scipy 
from io import StringIO
import statistics
import fileinput
import sys
import random
import time
import os


def fill_mat(matrix, sb):
	i = 0
	for line in sb:
		matrix[i] = line.strip().split()
		i = i + 1
	matrix[0].insert(0,'-')
	return matrix


def fill_s(s):
	
	seq = []
	boolean = 0
	
	for line in s:
		if (boolean != 0):
			i = 0
			while (i < len(line) and line[i] != '\n'):
				seq.append(line[i])
				i = i+1
		boolean = 1
	
	s.close()
	return seq


# find score in the substitution matrix
def find_score(v1,v2):
	
	col = matrix[0].index(v1)
	for i in range(1, 5):
		try:
			x = matrix[i].index(v2)
		except ValueError:
			continue;
		else:
			row = i
			break
	
	return matrix[row][col]


# init the 2D matrix of dim with 0s
def mat_init_zeros(dim):
	
	result = []
	for x in range(dim[0]):
		result.append([])
		for y in range(dim[1]):
			result[-1].append(0)
	
	return result


def calculate_currow(prev_row, n, seq1, seq2, prev_row_idx):
	curr_row = mat_init_zeros((1, n+1))
	curr_row[0][0] = prev_row[0][0] + gap_penalty
	for k in range(1, n+1):
		insert_horizontal = curr_row[0][k-1] + gap_penalty
		delete_vertical = prev_row[0][k] + gap_penalty
		match_diagonal = prev_row[0][k-1] + int(find_score(seq1[prev_row_idx], seq2[k-1]))
		curr_row[0][k] = max(insert_horizontal, delete_vertical, match_diagonal)
	return curr_row


def godown_nextrow(prev_row, curr_row, n, seq1, seq2, prev_row_idx):
	next_row = calculate_currow(curr_row, n, seq1, seq2, prev_row_idx+1)
	return curr_row,next_row


# DP align two seqs
def align(seq1, seq2):
	
	m = len(seq1)
	n = len(seq2)

	## DP memorization Score table of size (2, n+1) -> prev_row + curr_row
	prev_row_idx = 0
	prev_row = mat_init_zeros((1, n+1))
	
	#Base Case
	for k in range(0, n+1): # 0到n 是n+1个数
		prev_row[0][k] = k * gap_penalty
	
	curr_row = calculate_currow(prev_row, n, seq1, seq2, prev_row_idx)

	#Induction Rule
	align1 = '' + seq1[0] # alignment for seq1
	align2 = '' + seq2[0]
	i = 1  # current row idx position
	j = 1  # current col idx position
	# prev_row starts at row 1
	prev_row, curr_row = godown_nextrow(prev_row, curr_row, n, seq1, seq2, prev_row_idx)
	prev_row_idx = prev_row_idx + 1	

	while (i < m and j < n and prev_row_idx < m-1):
		#三种可能走法，从当前点，右/下/右下
		score_godown = curr_row[0][j]
		score_goright = prev_row[0][j + 1]
		score_gorightdown = curr_row[0][j + 1]
		best_movement = max(score_godown, score_goright, score_gorightdown)
		#insert_horizontal
		if (best_movement == score_goright):
			align1 = align1 + '-'
			align2 = align2 + seq2[j]
			j = j + 1
			continue;
		#delete_vertical
		elif (best_movement == score_godown):
			align1 = align1 + seq1[i]
			align2 = align2 + '-'
			i = i + 1
			prev_row, curr_row = godown_nextrow(prev_row, curr_row, n, seq1, seq2, prev_row_idx)
			prev_row_idx = prev_row_idx + 1	
			continue;
		#match_diagonal
		elif (best_movement == score_gorightdown):
			align1 = align1 + seq1[i]
			align2 = align2 + seq2[j]
			i = i + 1
			j = j + 1	
			prev_row, curr_row = godown_nextrow(prev_row, curr_row, n, seq1, seq2, prev_row_idx)
			prev_row_idx = prev_row_idx + 1	
			continue;				
	
	while (i < m and prev_row_idx < m-1):
		align1 = align1 + seq1[i]
		align2 = align2 + '-'
		i = i + 1
		prev_row, curr_row = godown_nextrow(prev_row, curr_row, n, seq1, seq2, prev_row_idx)
		prev_row_idx = prev_row_idx + 1	
	
	while (j < n):
		align1 = align1 + '-'
		align2 = align2 + seq2[j]
		j = j + 1


	best_score = curr_row[0][n]

	print ("The optimal alignment between given sequences has score " + str(best_score) + '.')
	print(align1)
	print(align2)
	
	filename = 'result.txt'
	f = open(filename,'w')
	f.write(align1)
	f.write("\n")
	f.write(align2)
	f.close()	


if __name__ == "__main__":
	
	s1 = open(sys.argv[1], 'r')
	s2 = open(sys.argv[2], 'r')
	sb = open(sys.argv[3], 'r')
	gap_penalty = int(sys.argv[4])
	
	# fill substitution matrix
	matrix = mat_init_zeros((5,5))
	matrix = fill_mat(matrix, sb)
		
	seq1 = fill_s(s1)	
	seq2 = fill_s(s2)
	
	align(seq1,seq2)
