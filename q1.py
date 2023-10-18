import blosum as bl
import numpy as np

mat = bl.BLOSUM(62)

s1 = "WKMDKSYWLFVREKKTDLCM"
n = len(s1)
s2 = "AIDDKSWAFVRECKTDQTW"
m = len(s2)
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
        _delete = scores[i-1][j] + indel_penalty
        _insert = scores[i][j-1] + indel_penalty

        # find the best-scoring choice
        # fill in score and backtrack matrices appropriately
        choice = max(match, _insert, _delete)
        if choice == match:
            backtrack[i][j] = 1
            scores[i][j] = match
        elif choice == _delete:
            backtrack[i][j] = -2
            scores[i][j] = _delete
        else:
            backtrack[i][j] = -1
            scores[i][j] = _insert
           
           
# recover the aligned strings via backtrack matrix  #
s1_aligned = ""
s2_aligned = ""
i = len(s1)
j = len(s2)

# start from sink of scoring matrix
# use backtrack arrows to return to source
#   1 = diagonal move
#  -2 = down move
#  -1 = right move
while (i > 0 or j > 0) and backtrack[i][j] != 0:
    # diagonal move
    if i > 0 and j > 0 and backtrack[i][j] == 1:
        s1_aligned = s1[i-1] + s1_aligned
        s2_aligned = s2[j-1] + s2_aligned
        i -= 1
        j -= 1
    # down move
    elif i > 0 and backtrack[i][j] == -2:
        s1_aligned = s1[i-1] + s1_aligned
        s2_aligned = "_" + s2_aligned
        i -= 1
    # right move
    elif j > 0 and backtrack[i][j] == -1:
        s1_aligned = "_" + s1_aligned
        s2_aligned = s2[j-1] + s2_aligned
        j -= 1

# perform any deletions at start if necessart
while i > 0:
    s1_aligned = s1[i-1] + s1_aligned
    s2_aligned = "_" + s2_aligned
    i -= 1

# perfomr any insertions at the start if necessary
while j > 0:
    s1_aligned = "_" + s1_aligned
    s2_aligned = s2[j-1] + s2_aligned
    j -= 1

print("score: " + str(scores[n][m]))
print(s1_aligned)
print(s2_aligned)