import sys
sys.setrecursionlimit(1500)


def backtrack_all_alignments(S, T, P, i, j, S_alg, T_alg, alignments):
    # Base case: if both sequences are aligned, store the result
    if i == 0 and j == 0:
        alignments.append(("".join(S_alg[::-1]), "".join(T_alg[::-1])))
        return

    # Recursively explore all possible directions (diagonal, left, top)
    if "d" in P[i][j]:  # Diagonal move
        S_alg.append(S[i - 1])
        T_alg.append(T[j - 1])
        backtrack_all_alignments(S, T, P, i - 1, j - 1, S_alg, T_alg, alignments)
        S_alg.pop()
        T_alg.pop()

    if "l" in P[i][j]:  # Left move
        S_alg.append("_")
        T_alg.append(T[j - 1])
        backtrack_all_alignments(S, T, P, i, j - 1, S_alg, T_alg, alignments)
        S_alg.pop()
        T_alg.pop()

    if "t" in P[i][j]:  # Top move
        S_alg.append(S[i - 1])
        T_alg.append("_")
        backtrack_all_alignments(S, T, P, i - 1, j, S_alg, T_alg, alignments)
        S_alg.pop()
        T_alg.pop()


# TODO: modify to include slippage once I validate that NW works
def align(S, T):
    m = len(S)
    n = len(T)

    # Initialize the scoring matrix
    A = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    P = [[[] for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(m + 1):
        P[i][0] = ["t"]

    for j in range(n + 1):
        P[0][j] = ["l"]

    P[0][0] = ["x"]

    # Fill the first row and first column with gap penalty
    for i in range(m + 1):
        A[i][0] = gap_penalty * i
    for j in range(n + 1):
        A[0][j] = gap_penalty * j

    # Fill the rest of the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if S[i - 1] == T[j - 1]:
                diag = A[i - 1][j - 1] + match_score
            else:
                diag = A[i - 1][j - 1] + mismatch_score

            top = A[i - 1][j] + gap_penalty
            left = A[i][j - 1] + gap_penalty

            max_val = max(diag, top, left)
            A[i][j] = max_val

            if diag == max_val:
                P[i][j] += ["d"]
            if top == max_val:
                P[i][j] += ["t"]
            if left == max_val:
                P[i][j] += ["l"]

    i, j = len(S), len(T)

    alignments = []
    backtrack_all_alignments(S, T, P, i, j, [], [], alignments)

    return (
        A[m][n],
        alignments,
    )  # Return all alignments and the final score (last matrix element)


def verify_alignment_score(S_aln, T_aln, calculated_score):
    assert len(S_aln) == len(T_aln)

    score = 0
    for i in range(len(S_aln)):
        if S_aln[i] == T_aln[i]:
            score += match_score
        elif S_aln[i] == "_" or T_aln[i] == "_":
            score += gap_penalty  # TODO: correct for slippage
        else:  # S_aln[i] != T_aln[i]
            score += mismatch_score

    if score != calculated_score:
        raise AssertionError("There's a mistake")


# TODO: change these or read from file (FASTA format)
# S = "GCAGGCAAGTGGGGCACCCGTATCCTTTCCAACTTACAAGGGTCCCCGTT"
# T = "GTGCGCCAGAGGAAGTCACTTTATATCCGCGCACGGTACTCCTTTTTCTA"

# S = "APPLE"
# T = "CAPE"

S = "TTCAGGCTGTTGTTGGCTTAGGGCTGGAAGCACAGAGTGGCTTGGCCTCAAGAGAATAGCTGGTTTCCCTAAGTTTACTTCTCTAAAACCCTGTGTTCACAAAGGCAGAGAGTCAGACCCTTCAATGGAAGGAGAGTGCTTGGGATCGATTATGTGACTTAAAGTCAGAATAGTCCTTGGGCAGTTCTCAAATGTTGGAGTGGAACATTGGGGAGGAAATTCTGAGGCAGGTATTAGAAATGAAAAGGAAACTTGAAACCTGGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGCAGATCACTGGAGGTCAGGAGTTCGAAACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAGAAATTAGCCGGTCATGGTGGTGGACACCTGTAATCCCAGCTACTCAGGTGGCTAAGGCAGGAGAATCACTTCAGCCCGGGAGGTGGAGGTTGCAGTGAGCCAAGATCATACCACGGCACTCCAGCCTGGGTGACAGTGAGACTGTGGCTCAAAAAAAAAAAAAAAAAAAGGAAAATGAAACTAGAAGAGATTTCTAAAAGTCTGAGATATATTTGCTAGATTTCTAAAGAATGTGTTCTAAAACAGCAGAAGATTTTCAAGAACCGGTTTCCAAAGACAGTCTTCTAATTCCTCATTAGTAATAAGTAAAATGTTTATTGTTGTAGCTCTGGTATATAATCCATTCCTCTTAAAATATAAGACCTCTGGCATGAATATTTCATATCTATAAAATGACAGATCCCACCAGGAAGGAAGCTGTTGCTTTCTTTGAGGTGATTTTTTTCCTTTGCTCCCTGTTGCTGAAACCATACAGCTTCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"
T = "TCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"

# S = "TGACGTGC"
# T = "TCGACGTCA"

match_score = 1
mismatch_score = -1
gap_penalty = -2
slippage_gap_penalty = None

score, alignments = align(S, T)
print("Alignment score:", score, "\n")

for i, alignment in enumerate(alignments):
    print("Alignment", i)
    print(alignment[0])
    print(alignment[1])
    verify_alignment_score(alignment[0], alignment[1], score)
    print()