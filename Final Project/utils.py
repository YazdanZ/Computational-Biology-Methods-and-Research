import random
from collections import defaultdict

import numpy as np

import PARAMS


def read_inputs(sequence_dir, confidence_dir):
    with open(sequence_dir, "r") as sequence_f:
        # strip removes any leading/trailing whitespace or newlines
        sequence = sequence_f.read().strip()

    with open(confidence_dir, "r") as conf_f:
        # Convert the numbers to a list of floats
        confidences = list(map(float, conf_f.read().split()))

    assert len(sequence) == len(
        confidences
    ), "Error: The length of sequence and confidence file do not match!"

    # probability array has a shape: (4, len(sequence))
    probability_array = np.zeros((4, len(sequence)))

    # same order as the handout
    nucleotide_to_row = {"A": 0, "C": 1, "G": 2, "T": 3}

    for i, (nucleotide, confidence) in enumerate(zip(sequence, confidences)):
        row_index = nucleotide_to_row[nucleotide]
        probability_array[row_index, i] = confidence

        # Distribute the remaining probability 1/3 among the other nucleotides
        incorrect_probability = (1 - confidence) / 3
        for j in range(4):
            if j != row_index:
                probability_array[j, i] = incorrect_probability

    return sequence, probability_array


def build_index(probability_array, w, threshold):
    n = probability_array.shape[1]
    index = defaultdict(list)

    # Fill in dictionary of seeds with confidence and indices
    for i in range(n - w + 1):
        wmer = ""
        confidence = 1

        for j in range(w):
            # Find the most likely nucleotide and its probability at each position
            position_probs = probability_array[:, i + j]
            max_row = np.argmax(position_probs)
            max_prob = position_probs[max_row]

            # Update the w-mer string and confidence score
            wmer += PARAMS.rows_to_nucleotides[max_row]
            confidence *= max_prob

        # Only index the w-mers with confidences > threshold
        if confidence >= threshold:
            index[wmer].append((i, confidence))

    return dict(index)


def add_substitutions(sequence, substitution_rate):
    num_substitutions = int(len(sequence) * substitution_rate / 100)

    for _ in range(num_substitutions):
        substitution_pos = random.randint(0, len(sequence) - 1)
        original_nucleotide = sequence[substitution_pos]
        possible_substitutions = ["A", "C", "G", "T"]
        possible_substitutions.remove(original_nucleotide)
        sequence[substitution_pos] = random.choice(possible_substitutions)

    return sequence


def add_indels(sequence, deletion_rate, insertion_rate):
    i = 0
    while i < len(sequence):
        if random.random() < deletion_rate / 100.0:
            sequence.pop(i)  # deletes element randomly
        elif random.random() < insertion_rate / 100.0:
            insertion = random.choice(["A", "C", "G", "T"])
            sequence.insert(i, insertion)
            i += 1
        i += 1
    return sequence


def gen_query_sequence(
    probability_array, query_length, substitution_rate, deletion_rate, insertion_rate, sanity=False
):
    # pick a random starting position
    start_pos = random.randint(0, probability_array.shape[1] - query_length)

    # generate query sequence based on probabilities at each position
    query_sequence = ""
    for i in range(query_length):
        position_probs = probability_array[:, start_pos + i]
        nucleotide = np.random.choice([0, 1, 2, 3], p=position_probs)
        if sanity:
            nucleotide = np.argmax(probability_array[:, start_pos + i])
        query_sequence += PARAMS.rows_to_nucleotides[nucleotide]

    query_sequence = list(
        query_sequence
    )  # following 2 functions called require lists instead od str
    query_sequence = add_substitutions(query_sequence, substitution_rate)
    query_sequence = add_indels(query_sequence, deletion_rate, insertion_rate)

    return "".join(query_sequence), start_pos


def score_alignment(query, reference, probability_array):
    score = 0.0  # initialize a score

    for i in range(len(query)):
        query_nucleotide = query[i]
        reference_nucleotide = reference[i]

        ref_idx = PARAMS.nucleotides_to_rows[reference_nucleotide]

        # Get the probability of the reference nucleotide at this position
        position_probs = probability_array[:, i]
        reference_confidence = position_probs[ref_idx]

        if query_nucleotide == reference_nucleotide:
            score += reference_confidence  # + score based on confidence of the match
        else:
            score -= reference_confidence  # - score based on reference confidence

    return score


# def find_probabilistic_matches(
#     query, probability_array, index, w, score_threshold, confidence_threshold
# ):
#     candidates = []

#     # Slide window over query sequence
#     for query_start in range(len(query) - w + 1):
#         query_wmer = query[query_start : query_start + w]

#         # Iterate over all w-mers in the index to allow partial matching
#         for ref_wmer, ref_entries in index.items():
#             # Calculate the probabilistic match score between the query and the reference w-mer
#             for ref_start, ref_confidence in ref_entries:
#                 # Skip low-confidence matches
#                 if ref_confidence < confidence_threshold:
#                     continue

#                 ref_probabilities = probability_array[:, ref_start : ref_start + w]

#                 # Compute the probabilistic match score
#                 raw_score = score_alignment(query_wmer, ref_wmer, ref_probabilities)

#                 # Weight the score by the reference confidence
#                 weighted_score = raw_score * ref_confidence

#                 # Only add to candidates if igger than the score threshold
#                 if weighted_score >= score_threshold:
#                     candidates.append((ref_start, weighted_score))

#     return candidates


def find_probabilistic_matches(
    query_sequence, probability_array, index, w, score_threshold, confidence_threshold
):
    candidates = []

    for query_start in range(len(query_sequence) - w + 1):
        query_wmer = query_sequence[query_start : query_start + w]

        # Check if the query w-mer exists in the index
        if query_wmer in index:
            for ref_start, ref_confidence in index[query_wmer]:
                # Skip low-confidence matches early
                if ref_confidence < confidence_threshold:
                    continue

                ref_wmer = "".join(
                    PARAMS.rows_to_nucleotides[
                        np.argmax(probability_array[:, ref_start + i])
                    ]
                    for i in range(w)
                )

                # Compute the probabilistic match score
                raw_score = score_alignment(
                    query_wmer,
                    ref_wmer,
                    probability_array[:, ref_start : ref_start + w],
                )

                # Weight the score by the reference w-mer confidence
                weighted_score = raw_score * ref_confidence

                if weighted_score >= score_threshold:
                    candidates.append((ref_start, weighted_score))

    return candidates


def ungapped_extension_probabilistic(
    query, probability_array, candidate_start, match_length, gap_penalty
):
    query_len = len(query)
    reference_len = probability_array.shape[1]

    alg_score = 0.0  # Init score

    # left extension
    left_extension = candidate_start
    for i in range(match_length):
        query_pos = match_length - 1 - i
        ref_pos = candidate_start - 1 - i

        if query_pos < 0 or ref_pos < 0:
            break  # Out of bounds

        ref_probs = probability_array[:, ref_pos]
        query_nucleotide = query[query_pos]

        match_prob = ref_probs[PARAMS.nucleotides_to_rows[query_nucleotide]]
        if query_nucleotide == PARAMS.rows_to_nucleotides[np.argmax(ref_probs)]:
            alg_score += match_prob  # reward reference confidence
        else:
            mismatch_prob = 1 - match_prob
            alg_score -= mismatch_prob

    # right extension
    right_extension = candidate_start + match_length
    for i in range(match_length, query_len):
        query_pos = i
        ref_pos = candidate_start + i

        if query_pos >= query_len or ref_pos >= reference_len:
            break  # Out of bounds

        ref_probs = probability_array[:, ref_pos]
        query_nucleotide = query[query_pos]

        match_prob = ref_probs[PARAMS.nucleotides_to_rows[query_nucleotide]]
        if query_nucleotide == PARAMS.rows_to_nucleotides[np.argmax(ref_probs)]:
            alg_score += match_prob
        else:
            mismatch_prob = 1 - match_prob
            alg_score -= mismatch_prob

    if left_extension > 0:
        alg_score -= gap_penalty * left_extension

    if right_extension < reference_len:
        alg_score -= gap_penalty * (query_len - right_extension)

    return alg_score


def gapped_extension(query, probability_array, candidate_start, gap_penalty):
    query_len = len(query)
    reference_len = probability_array.shape[1]

    # Limit the alignment range to avoid out-of-bounds access
    max_j = min(reference_len - candidate_start, query_len)

    # Initialize scoring matrix and traceback matrix
    dp = np.zeros((query_len + 1, max_j + 1))
    traceback = np.zeros((query_len + 1, max_j + 1), dtype=int)

    DIAG, UP, LEFT = 1, 2, 3

    for i in range(1, query_len + 1):
        dp[i, 0] = dp[i - 1, 0] - gap_penalty
        traceback[i, 0] = UP
    for j in range(1, max_j + 1):
        dp[0, j] = dp[0, j - 1] - gap_penalty
        traceback[0, j] = LEFT

    for i in range(1, query_len + 1):
        for j in range(1, max_j + 1):
            query_nucleotide = query[i - 1]
            ref_pos = candidate_start + j - 1

            # for ref_pos to not exceed probability_array bounds
            if ref_pos >= reference_len:
                continue

            match_prob = probability_array[
                PARAMS.nucleotides_to_rows[query_nucleotide], ref_pos
            ]
            mismatch_prob = 1 - match_prob

            diag_score = dp[i - 1, j - 1] + (
                match_prob
                if query_nucleotide
                == PARAMS.rows_to_nucleotides[np.argmax(probability_array[:, ref_pos])]
                else -mismatch_prob
            )
            up_score = dp[i - 1, j] - gap_penalty
            left_score = dp[i, j - 1] - gap_penalty

            dp[i, j] = max(diag_score, up_score, left_score)

            # Traceback
            if dp[i, j] == diag_score:
                traceback[i, j] = DIAG
            elif dp[i, j] == up_score:
                traceback[i, j] = UP
            else:
                traceback[i, j] = LEFT

    # TRACEBACK
    aligned_query = []
    aligned_reference = []

    i, j = query_len, max_j
    while i > 0 or j > 0:
        if traceback[i, j] == DIAG:
            aligned_query.append(query[i - 1])
            aligned_reference.append(
                PARAMS.rows_to_nucleotides[
                    np.argmax(probability_array[:, candidate_start + j - 1])
                ]
            )
            i -= 1
            j -= 1
        elif traceback[i, j] == UP:
            aligned_query.append(query[i - 1])
            aligned_reference.append("-")
            i -= 1
        elif traceback[i, j] == LEFT:
            aligned_query.append("-")
            aligned_reference.append(
                PARAMS.rows_to_nucleotides[
                    np.argmax(probability_array[:, candidate_start + j - 1])
                ]
            )
            j -= 1

    aligned_query.reverse()
    aligned_reference.reverse()

    return dp[query_len, max_j], "".join(aligned_query), "".join(aligned_reference)

