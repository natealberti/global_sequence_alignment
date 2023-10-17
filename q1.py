import blosum as bl
import numpy as np

mat = bl.BLOSUM(62)

s1 = "WKMDKSYWLFVREKKTDLCM"
n = len(s1)
s2 = "AIDDKSWAFVRECKTDQTW"
m = len(s2)
indel_penalty = 5

# initialize n x m matrix of zeros
scores = np.zeros((n, m))
backtrack = np.zeros((n, m))

for i in range(1, n):
    for j in range(1, m):
        match = scores[i-1][j-1] + mat[s1[i-1]][s2[j-1]]
        _delete = scores[i-1][j] - indel_penalty
        _insert = scores[i][j-1] - indel_penalty
        e = max(match, _delete, _insert)
        if e == _delete:
            backtrack[i][j] = -2
        elif e == _insert:
            backtrack[i][j] = -1
        else:
            backtrack[i][j] = 1
        scores[i][j] = max(
            scores[i-1][j-1] + mat[s1[i-1]][s2[j-1]],   # match/mismatch
            scores[i-1][j] - indel_penalty,             # delete
            scores[i][j-1] - indel_penalty              # insert
        )

def recoverLCS(backtrack, s1, s2):
    s1_aligned = ""
    s2_aligned = ""
    i = len(s1)-1
    j = len(s2)-1
    s1_accumulator = 0
    s2_accumulator = 0
    while i > 0 and j > 0:
        if backtrack[i][j] == -2:
            s1_aligned += s1[s1_accumulator]
            s1_accumulator += 1
            s2_aligned += "_"
            i -= 1
        elif backtrack[i][j] == -1:
            s2_aligned += s2[s2_accumulator]
            s2_accumulator += 1
            s1_aligned += "_"
            j -= 1
        else:
            s1_aligned += s1[s1_accumulator]
            s1_accumulator += 1
            s2_aligned += s2[s2_accumulator]
            s2_accumulator += 1
            i -= 1
            j -= 1
    return [s1_aligned, s2_aligned]

print("score: " + str(scores[n-1][m-1]))
s1_aligned, s2_aligned = recoverLCS(backtrack, s1, s2)
print(s1_aligned)
print(s2_aligned)