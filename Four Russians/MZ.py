import numpy as np 
import itertools,math

t = 2


def get_binary_differences(s):
    s = tuple(int(m-n) for n,m in zip(s,s[1:]))
    return(s)


def preprocessing():
	# precompute all the possible inputs in terms of variables{A, B, 0, C,E} 
	                                       # choices = [-1, 0, 1] 
	#[(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
	A_possibility = list(itertools.product(range(-1,2), repeat = t))
	B_possibility = list(itertools.product(range(-1,2), repeat = t))
	#['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
	C_possibility = ["".join(x) for x in itertools.product("ACGT", repeat = t)]
	E_possibility = ["".join(x) for x in itertools.product("ACGT", repeat = t)]

	LUP = {}
	for e1 in A_possibility:
		for e2 in B_possibility:
			for e3 in C_possibility:
				for e4 in E_possibility:
					A_, B_, K_ = block_compute(e1, e2, e3, e4)
					LUP[(e1,e2,e3,e4)] = (A_,B_,K_)

	return LUP[((-1,-1),(0,1),'AG','CC')]


def block_compute(A, B, C, E):
	# memorization table
	DP = np.zeros((t+1, t+1))

	# base case
	
	for i in range(1, t+1):
		DP[i][0] = A[i-1]
		DP[0][i] = B[i-1]

	for row in range(1,t+1):
		for col in range(1,t+1):
			diag = 0
			if C[row-1] == E[col-1]:
				diag = 0
			else:
				diag = 1

			DP[row][col] = min(min(DP[row-1][col]+1 ,DP[row][col-1]+1), DP[row-1][col-1] + diag)

	# (A_,B_,K_)
	A_ = get_binary_differences(DP[:,t])
	B_ = get_binary_differences(DP[t,:])
	K_ = DP[t][t]

	return A_,B_,K_

if __name__ == '__main__':
    print(preprocessing())









