import blosum as bl
import numpy as np
import sys

mat = bl.BLOSUM(62)

# global alignment function
def local_alignment(s1, s2):
    n, m = len(s1), len(s2)
    indel_penalty = -11
    epsilon_penalty = -1
    
    # initialize scoring and backtrack matrices
    middle = np.zeros((n+1, m+1))
    upper = np.zeros((n+1, m+1))
    lower = np.zeros((n+1, m+1))

    # fill scoring matrix via dynamic programming
    for i in range(1, n+1):
        for j in range(1, m+1):
            # find the score of each cell with match, insert, delete
            # lower matrix recurrence
            lower[i][j] = max(lower[i-1][j] + epsilon_penalty, middle[i-1][j] + indel_penalty)
            
            # upper matrix recurrence
            upper[i][j] = max(upper[i][j-1] + epsilon_penalty, middle[i][j-1] + indel_penalty)
            
            # middle matrix recurrence
            middle[i][j] = max(lower[i][j], middle[i-1][j-1] + mat[s1[i-1]][s2[j-1]], upper[i][j], 0)
            
    # recover the aligned strings via backtrack matrix  #
    s1_aligned, s2_aligned = "", ""

    # find largest value in middle matrix and its indices
    i_max, j_max = np.unravel_index(np.argmax(middle), middle.shape)
    
    # start from the best-scoring node
    i, j = i_max, j_max
    while (i > 0 or j > 0) and middle[i][j] > 0:
        # diagonal move
        if i > 0 and j > 0 and middle[i][j] != upper[i][j] and middle[i][j] != lower[i][j]:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            i -= 1
            j -= 1
        # right move
        elif j > 0 and middle[i][j] == upper[i][j]:
            s1_aligned = "" + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            j -= 1  
        # down move
        elif i > 0 and middle[i][j] == lower[i][j]:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = "" + s2_aligned
            i -= 1
    
    return [s1_aligned, s2_aligned, middle[i_max][j_max]]

# find the file to open
try:
    file_name = sys.argv[1]
except:
    file_name = "test_1.txt"
    
# open the file
file = open(file_name, 'r')
file.readline()
s1, s2 = "", ""
next_seq = False
for line in file:
    line = line.replace("\n", "")
    if line != ">Sequence 2" and not next_seq:
        s1 += line
    elif line != ">Sequence 2" and next_seq:
        s2 += line
    else:
        next_seq = True
    
# perform local alignment
s1_aligned, s2_aligned, score = local_alignment(s1, s2)
w = open("output_q2_Nathaniel_Alberti.txt", "w")
w.writelines(str(int(score)) + "\n")
w.write(s1_aligned + "\n")
w.write(s2_aligned)