import sys

sys.setrecursionlimit(3000)


def backtrack_all_alignments(S, T, P, i, j, S_alg, T_alg, alignments, find_one=False):
    # Base case: if both sequences are aligned, store the result
    if i == 0 and j == 0:
        alignments.append(("".join(S_alg[::-1]), "".join(T_alg[::-1])))
        # print("".join(S_alg[::-1]), "".join(T_alg[::-1]))
        return

    # skip recursion if we have an alignment and we only want to find one
    if find_one and alignments:
        return

    # Recursively explore all possible directions (diagonal, left, top)
    if "d" in P[i][j]:  # Diagonal move
        S_alg.append(S[i - 1])
        T_alg.append(T[j - 1])
        backtrack_all_alignments(
            S, T, P, i - 1, j - 1, S_alg, T_alg, alignments, find_one
        )
        S_alg.pop()
        T_alg.pop()

    if find_one and alignments:
        return

    if "l" in P[i][j]:  # Left move
        S_alg.append("_")
        T_alg.append(T[j - 1])
        backtrack_all_alignments(S, T, P, i, j - 1, S_alg, T_alg, alignments, find_one)
        S_alg.pop()
        T_alg.pop()

    if find_one and alignments:
        return

    if "t" in P[i][j]:  # Top move
        S_alg.append(S[i - 1])
        T_alg.append("_")
        backtrack_all_alignments(S, T, P, i - 1, j, S_alg, T_alg, alignments, find_one)
        S_alg.pop()
        T_alg.pop()


def align(S, T, find_one=False):
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

            if S[i - 1] == S[i - 2]:
                top = A[i - 1][j] + slippage_gap_penalty
            else:
                top = A[i - 1][j] + gap_penalty

            if T[j - 1] == T[j - 2]:
                left = A[i][j - 1] + slippage_gap_penalty
            else:
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

    # for troubleshooting
    # for i in range(len(A)):
    #     print(A[i])
    # for i in range(len(P)):
    #     print(P[i])

    alignments = []
    backtrack_all_alignments(S, T, P, i, j, [], [], alignments, find_one)

    return (
        A[m][n],
        alignments,
    )  # Return all alignments and the final score (last matrix element)


def verify_alignment_score(S_aln, T_aln, calculated_score):
    assert len(S_aln) == len(T_aln)

    score = 0
    for i in range(len(S_aln)):
        if i != 0:  # so that it doesn't compare element 0 to last element when i=0
            if S_aln[i] == T_aln[i]:
                score += match_score
            elif S_aln[i] == "_" and T_aln[i] == T_aln[i - 1]:
                score += slippage_gap_penalty
            elif S_aln[i] == "_" and T_aln[i] != T_aln[i - 1]:
                score += gap_penalty
            elif T_aln[i] == "_" and S_aln[i] == S_aln[i - 1]:
                score += slippage_gap_penalty
            elif T_aln[i] == "_" and S_aln[i] != S_aln[i - 1]:
                score += gap_penalty
            else:  # S_aln[i] != T_aln[i]
                score += mismatch_score
        else:
            if S_aln[i] == T_aln[i]:
                score += match_score
            elif S_aln[i] == "_" or T_aln[i] == "_":
                score += gap_penalty
            else:
                score += mismatch_score

    if score != calculated_score:
        raise AssertionError(
            f"Expected score:{score} doesn't match alignment score {calculated_score}"
        )


def get_sequence_from_fasta(file_path, name):
    with open(file_path, 'r') as file:
        found_name = False
        for line in file:
            line = line.strip()  # Clean up whitespace
            if found_name:
                return line  # Return the line directly below the header
            if line == f">{name}":
                found_name = True
    return "Sequence not found."



# TODO: change these or read from file (FASTA format)
# S = "GCAGGCAAGTGGGGCACCCGTATCCTTTCCAACTTACAAGGGTCCCCGTT"
# T = "GTGCGCCAGAGGAAGTCACTTTATATCCGCGCACGGTACTCCTTTTTCTA"

# S = "APPLE"
# T = "CAPE"

# S = "TTCAGGCTGTTGTTGGCTTAGGGCTGGAAGCACAGAGTGGCTTGGCCTCAAGAGAATAGCTGGTTTCCCTAAGTTTACTTCTCTAAAACCCTGTGTTCACAAAGGCAGAGAGTCAGACCCTTCAATGGAAGGAGAGTGCTTGGGATCGATTATGTGACTTAAAGTCAGAATAGTCCTTGGGCAGTTCTCAAATGTTGGAGTGGAACATTGGGGAGGAAATTCTGAGGCAGGTATTAGAAATGAAAAGGAAACTTGAAACCTGGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGCAGATCACTGGAGGTCAGGAGTTCGAAACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAGAAATTAGCCGGTCATGGTGGTGGACACCTGTAATCCCAGCTACTCAGGTGGCTAAGGCAGGAGAATCACTTCAGCCCGGGAGGTGGAGGTTGCAGTGAGCCAAGATCATACCACGGCACTCCAGCCTGGGTGACAGTGAGACTGTGGCTCAAAAAAAAAAAAAAAAAAAGGAAAATGAAACTAGAAGAGATTTCTAAAAGTCTGAGATATATTTGCTAGATTTCTAAAGAATGTGTTCTAAAACAGCAGAAGATTTTCAAGAACCGGTTTCCAAAGACAGTCTTCTAATTCCTCATTAGTAATAAGTAAAATGTTTATTGTTGTAGCTCTGGTATATAATCCATTCCTCTTAAAATATAAGACCTCTGGCATGAATATTTCATATCTATAAAATGACAGATCCCACCAGGAAGGAAGCTGTTGCTTTCTTTGAGGTGATTTTTTTCCTTTGCTCCCTGTTGCTGAAACCATACAGCTTCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"
# T = "TCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"

# S = "TGACGTGC"
# T = "TCGACGTCA"

# S = "AAAG"
# T = "AT"

# S = "ATTAGT"
# T = "ATAAG"

# S = "ATGA"
# T = "ACCGGGA"

# S = "AT"
# T = "AACG"

# S = "AGAAGGGGAA"
# T = "GAAA"

match_score = 1
mismatch_score = -1
gap_penalty = -2
slippage_gap_penalty = -1
find_all_alignments = False

file_path = 'hw1_brca1_3utr_small.fa'
S = get_sequence_from_fasta(file_path, "human_BRCA1_3UTR")
T = get_sequence_from_fasta(file_path, "mouse_BRCA1_3UTR")


score, alignments = align(S, T, find_one=not find_all_alignments)

with open('hw1_brca1_3utr_small_output.txt', 'w') as f:
    f.write(f"Alignment score: {score} \n\n")

    for i, alignment in enumerate(alignments):
        # f.write(f"Alignment {i + 1} \n")
        f.write(f"{alignment[0]} \n\n")
        f.write(f"{alignment[1]} \n\n")
        verify_alignment_score(alignment[0], alignment[1], score)

# print("Alignment score:", score, "\n")

# for i, alignment in enumerate(alignments):
#     print("Alignment", i + 1)
#     print(alignment[0])
#     print(alignment[1])
#     verify_alignment_score(alignment[0], alignment[1], score)
#     print()
