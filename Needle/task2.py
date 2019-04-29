#create_data_and_align  program should take as input a single integer L
#It should then create two “related” DNA sequences

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
import math

def mutate_seq(S_i):
	for i in range(0, int(L/10)):
		#mutate index, randomly choose a position in the sequence
		idx = random.randrange(0, len(S_i), 1)
		b = random.randrange(0, 2, 1)
		
		if (b == 1):
			#mutate it to something else with probability ½
			choices = list(mole_type)
			#exclude the possibility that A mutate to A itself
			choices.pop(mole_type.index(S_i[idx]))
			b = random.randrange(0, 3, 1)
			S_i[idx] = choices[b]
		
		elif (b == 0):
			#delete that position with probability ½
			S_i.pop(idx)
	
	return ''.join(S_i)


if __name__ == "__main__":
	
	mole_type = ['A','T','C','G']
	#use the current time to seed your pseudo-random number generator
	random.seed(time.time())
	
	#Generate S1
	S1 = []
	L = int(sys.argv[1])
	for i in range(0, L):
		c = random.randrange(0, 4, 1)
		S1.append(mole_type[c])
	
	#Generate S2 and S3
	for i in range(2,4):
		seq = list(S1)
		new_seq = mutate_seq(seq)
		#write/store to local file
		filename = 'S' + str(i) + '.txt'
		f = open(filename, 'w')
		f.write(">" + filename + "\n")
		f.write(new_seq)
		f.close()

# The “create_data_and_align” program should then call your alignment program 
# from the first task, and thus produce a new text file result.txt that has the optimal 
# alignment of sequences S2 and S3. Call your alignment program using the subs.txt 
# scoring matrix mentioned above and gap penalty = -500.
	cmd = "python task1.py S2.txt S3.txt subs.txt -500"
	os.system(cmd)
