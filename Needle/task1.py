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


# init the matrix of dim with 0s
def mat_init_zeros(dim):
	
	result = []
	for x in range(dim[0]):
		result.append([])
		for y in range(dim[1]):
			result[-1].append(0)
	
	return result


# DP align two seqs
def align(seq1, seq2):
	
	m = len(seq1)
	n = len(seq2)

	## DP memorization Score table
	score = mat_init_zeros((m+1, n+1))
	# m row seq1 * n col seq2
	
	#Base Case
	#score[0][0] is score[*][*] = 0
	for i in range(0, m+1):
		score[i][0] = i * gap_penalty
	for i in range(0, n+1):
		score[0][i] = i * gap_penalty
	
	#Induction Rule
	for i in range(1, m+1): #从上到下
		for j in range(1, n+1): #从左到右
			insert_horizontal = score[i][j-1] + gap_penalty #insert a '-' to seq1
			delete_vertical = score[i-1][j] + gap_penalty	#insert a '-' to seq2
			match_diagonal = score[i-1][j-1] + int(find_score(seq1[i-1], seq2[j-1]))
			score[i][j] = max(insert_horizontal, delete_vertical, match_diagonal)
	
	## Traceback reveal the path
	i = m
	j = n
	align1 = ''
	align2 = ''

	while (i > 0 and j > 0):
		cur = score[i][j]
		zuo = score[i-1][j]
		shang = score[i][j-1]
		zuo_shang = score[i-1][j-1]
				
		if (cur == zuo + gap_penalty):
			align1 = align1 + seq1[i-1]
			align2 = align2 + '-'
			i = i-1

		elif (cur == shang + gap_penalty):
			align1 = align1 + '-'
			align2 = align2 + seq2[j-1]
			j = j-1
		
		elif (cur == zuo_shang + int(find_score(seq1[i-1], seq2[j-1]))):
			align1 = align1 + seq1[i-1] # is the last elem in seq1, it is also score[i][j] represents
			align2 = align2 + seq2[j-1]
			i = i-1
			j = j-1		

	
	while (i > 0):
		align1 = align1 + seq1[i-1]
		align2 = align2 + '-'
		i = i-1
	
	while (j > 0):
		align1 = align1 + '-'
		align2 = align2 + seq1[j-1]
		j = j-1
	
	print ("The optimal alignment between given sequences has score " + str(score[m][n]) + '.')
	print(align1[::-1])
	print(align2[::-1])
	
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

	align("AGGCATGGC", "TAGCTATCA")
