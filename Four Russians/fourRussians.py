import numpy as np 
import itertools,math
import pickle,os
import time,sys

t = 3
backtracking = {}

def block_backtrack(s,t,i,j,T):
    #Generates aligned strings based on different starting points
    s_align, t_align = "", ""
    while i > 0 and j > 0:
        if s[i] == t[j] and T[i, j] == T[i-1, j - 1]: #match
            s_align = s[i] + s_align
            t_align = t[j] + t_align
            j -= 1
            i -= 1
        elif T[i, j] == T[i, j - 1] + 1:
            s_align = "_" + s_align
            t_align = t[j] + t_align
            j -= 1
        elif T[i, j] == T[i - 1, j] + 1:
            s_align = s[i] + s_align
            t_align = "_" + t_align
            i -= 1
        else:
            s_align = s[i] + s_align
            t_align = t[j] + t_align
            j -= 1
            i -= 1
    return((i,j),(str(s_align),str(t_align)))


def get_differences(s):
    s = tuple(int(m-n) for n,m in zip(s,s[1:]))
    return(s)


def preprocessing():
	# precompute all the possible inputs in terms of variables{A, B, 0, C,E} 
	                                       # choices = [-1, 0, 1] 
	#[(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
    A_possibility = list(itertools.product(range(-1,2), repeat = t-1))
    B_possibility = list(itertools.product(range(-1,2), repeat = t-1))
	#['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    C_possibility = ["".join(x) for x in itertools.product("ACGT", repeat = t)]
    E_possibility = ["".join(x) for x in itertools.product("ACGT", repeat = t)]

    C_possibility_emptyStr = ["_"+"".join(x) for x in itertools.product("ACGT", repeat = t-1)]
    E_possibility_emptyStr = ["_"+"".join(x) for x in itertools.product("ACGT", repeat = t-1)]

    C_possibility += C_possibility_emptyStr
    E_possibility += E_possibility_emptyStr

    LUP = {}
    for e1 in A_possibility:
	    for e2 in B_possibility:
			for e3 in C_possibility:
				for e4 in E_possibility:
					(A_,B_,back_tracking_ht) = block_compute(e1, e2, e3, e4)
					LUP[(e1,e2,e3,e4)] = (A_,B_,back_tracking_ht)

    return LUP



def block_compute(A, B, C, E):
    # memorization table
    DP = np.zeros((t, t))
    DP[0][0] = 0
    for i in range (1, t):
        DP[0, i] = B[i-1]+DP[0,i-1]
        DP[i, 0] = A[i-1]+DP[i-1,0]
 

    for i in range(1, t):
        for j in range(1, t):
            diag = 0
            if C[i] != E[j]:
                diag = 1
            
            DP[i][j] = min(min(DP[i,j-1]+1, DP[i-1,j]+1), DP[i-1,j-1]+diag)
    tmp = DP[:,t-1]

    A_ = tuple(int(m-n) for n,m in zip(tmp,tmp[1:]))
    tmp = DP[t-1,:]
    B_ = tuple(int(m-n) for n,m in zip(tmp,tmp[1:]))
    if 2 in A_ or 2 in B_:
        print("")
    

    back_tracking_ht = {}
                    
    #generate all possible alignments and endpoints depending on entry point
    for i in range(1,t):
        j=t-1
        (coords),(alignments)=block_backtrack(C,E,i,j,DP)
        back_tracking_ht[i,j] = (coords,alignments)
    for j in range(1,t):
        i =t-1
        (coords), (alignments) = block_backtrack(C, E, i, j, DP)
        back_tracking_ht[i,j] = (coords,alignments)
    return (A_, B_, back_tracking_ht)
            

        
def main(str1, str2):
    #if lookup table exists, just read it; else, create it and store
    exists = os.path.isfile('dict')
    LUP = {}
    if exists:
        with open("dict", "rb") as file:
            LUP=pickle.load(file)
    else:
        LUP = preprocessing()
        with open("dict", "wb") as file:
            pickle.dump(LUP, file)
    
    
    #initialize dp array
    n = len(str1)

    DP_final = np.zeros((n+1, n+1))
    DP_final[0] = [i for i in range(n+1)]
    DP_final[:,0]=[i for i in range(n+1)]

    
    str1_new = "_"+str1
    str2_new = "_"+str2

    start_time = time.time()
    #officially start main loop: fill dp table block by block
    for i in range(0, n, t-1):
        for j in range(0, n, t-1):
            A = DP_final[i:i+t, j]
            B = DP_final[i,j:j+t]
            #k = DP_final[i][j]
            
            e1 = tuple(j-i for i, j in zip(A[:-1], A[1:]))
            e2 = tuple(j-i for i, j in zip(B[:-1], B[1:]))
            e3 = str1_new[i:i+t]
            e4 = str2_new[j:j+t]

            
            (A_,B_,bt) = LUP[(e1,e2,e3,e4)] 

            for index in range(0, t-1):
                DP_final[i+t-1, j+1+index] = B_[index] + DP_final[i+t-1,j+index]
                DP_final[i+index+1, j+t-1] = A_[index] + DP_final[i+index, j+t-1]

    print("--- %s seconds ---" % (time.time() - start_time))

    #print DP_final
    final_s1, final_s2 = traceback(str1, str2, DP_final, LUP)
    print "score: " + str(DP_final[len(DP_final)-1, len(DP_final)-1])
    print "s1: "+final_s1
    print "s2: "+final_s2




def traceback(str1, str2, DP_final, LUP):
    i = len(DP_final)-1
    j = len(DP_final)-1

    str1_new = "_"+str1
    str2_new = "_"+str2

    index_map = {}
    for index in range(0, len(DP_final)-1, t-1):
        this_range = tuple(range(index,index+t))
        for each in range(index, index+t):
            if each in index_map:
                index_map[each] = (index_map[each], this_range)
            else:
                index_map[each] = this_range
    
    #print index_map
    #{0: (0, 1, 2), 
    # 1: (0, 1, 2), 
    # 2: ((0, 1, 2), (2, 3, 4)), 
    # 3: (2, 3, 4), 
    # 4: ((2, 3, 4), (4, 5, 6)), 
    # 5: (4, 5, 6), 
    # 6: (4, 5, 6)}

    final_s1 = ""
    final_s2 = ""
    
    while(i >= 0 and j >= 0):
        if i == 0 and j == 0:
            break

        if i in range(0, t) and j == 0:
            final_s1 = str1_new[i] + final_s1
            final_s2 = "_" + final_s2
            i = i -1
            continue
        
        elif i == 0 and j in range(0, t):
            final_s1 = "_" + final_s1
            final_s2 = str2_new[j] + final_s2
            j = j -1
            continue

        which_row = index_map[i]
        which_col = index_map[j]

        if not isinstance(which_row[0],int) and not isinstance(which_col[0],int):
            for each_r in which_row:
                ti = each_r.index(i)
                for each_c in which_col:
                    tj = each_c.index(j)
                    if ti != 0 and tj != 0: #not first row or column
                        which_row = each_r
                        which_col = each_c
                        break
                else:
                    continue
                break

        elif not isinstance(which_row[0],int):
            for each_r in which_row:
                ti = each_r.index(i)
                if ti != 0:
                    which_row = each_r
                    break

        elif not isinstance(which_col[0],int):
            for each_c in which_col:
                tj = each_c.index(j)
                if tj != 0:
                    which_col = each_c
                    break
    

        A = DP_final[which_row[0]:which_row[0]+t, which_col[0]]
        B = DP_final[which_row[0], which_col[0]:which_col[0]+t]
      
        e1 = get_differences(A)
        e2 = get_differences(B)
        e3 = str1_new[which_row[0]:which_row[0]+t]
        e4 = str2_new[which_col[0]:which_col[0]+t]

        (A_,B_,back_tracking_ht) = LUP[(e1,e2,e3,e4)] 
        #back_tracking_ht[i,j] = (coords,alignments) -> ((i,j),(s_align,t_align))
        block_i = which_row.index(i)
        block_j = which_col.index(j)

        ((i,j),(s_align,t_align)) = back_tracking_ht[block_i, block_j]
        final_s1 = s_align + final_s1
        final_s2 = t_align + final_s2

        new_i = which_row[i]
        new_j = which_col[j]

        i = new_i
        j = new_j


    return final_s1, final_s2


if __name__ == '__main__':
  
    #traceback("", "", np.zeros((7,7)), None)
    #main("AGTCGCAGGTTACCGT", "GCTAGGATCATGCAGT")
    #main("AGCGAA", "GCTAGA")
    #print(block_compute((1,1), (1, 1), '*GC', '*AG'))
    if sys.argv[1] and sys.argv[2]:
        x=""
        y=""
        s1name = sys.argv[1]
        s2name = sys.argv[2]

        s1File = open(s1name, "r")
        for line in s1File:
            if line[0] != '>':
			    x += line


        s2File = open(s2name, "r")
        for line in s2File:
            if line[0] != '>':
                y += line

        main(x,y)

    


