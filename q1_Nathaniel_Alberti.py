import blosum as bl
import numpy as np
import sys

mat = bl.BLOSUM(62)

# global alignment function
def global_alignment(s1, s2):
    n, m = len(s1), len(s2)
    indel_penalty = -5
    
    # initialize scoring and backtrack matrices
    scores = np.zeros((n+1, m+1))
    backtrack = np.zeros((n+1, m+1))

    # fill in first row and column of scoring matrix with indel penalties
    for i in range(1, n+1):
        scores[i][0] = indel_penalty*i
    for j in range(1, m+1):
        scores[0][j] = indel_penalty*j

    # fill scoring matrix via dynamic programming
    for i in range(1, n+1):
        for j in range(1, m+1):
            # find the score of each cell with match, insert, delete
            match = scores[i-1][j-1] + mat[s1[i-1]][s2[j-1]]
            _insert = scores[i-1][j] + indel_penalty
            _delete = scores[i][j-1] + indel_penalty

            # find the best-scoring choice
            # fill in score and backtrack matrices appropriately
            choice = max(match, _insert, _delete)
            if choice == match:
                backtrack[i][j] = 1
                scores[i][j] = match
            elif choice == _delete:
                backtrack[i][j] = -1
                scores[i][j] = _delete
            else:
                backtrack[i][j] = -2
                scores[i][j] = _insert
            
    # recover the aligned strings via backtrack matrix  #
    s1_aligned, s2_aligned = "", ""
    i, j = len(s1), len(s2)

    # start from sink of scoring matrix
    # use backtrack arrows to return to source
    #   1 = diagonal move
    #  -2 = down move
    #  -1 = right move
    while (i > 0 or j > 0) and backtrack[i][j] != 0:
        # diagonal move / match / mismatch
        if i > 0 and j > 0 and backtrack[i][j] == 1:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            i -= 1
            j -= 1
        # right move / deletion
        elif j > 0 and backtrack[i][j] == -1:
            s1_aligned = "-" + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            j -= 1
        # down move / insertion
        elif i > 0 and backtrack[i][j] == -2:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = "-" + s2_aligned
            i -= 1

    # perform any insertions at start if necessart
    while i > 0:
        s1_aligned = s1[i-1] + s1_aligned
        s2_aligned = "-" + s2_aligned
        i -= 1

    # perfomr any deletions at the start if necessary
    while j > 0:
        s1_aligned = "-" + s1_aligned
        s2_aligned = s2[j-1] + s2_aligned
        j -= 1
    
    return [s1_aligned, s2_aligned, scores[n][m]]

#  open file
# parse variables
file = open(sys.argv[1], 'r')
_ = file.readline()
s1 = file.readline().replace("\n", "")
_ = file.readline()
s2 = file.readline().replace("\n", "")

# perform global alignment
s1_aligned, s2_aligned, score = global_alignment(s1, s2)
print(str(int(score)))
print(s1_aligned)
print(s2_aligned)