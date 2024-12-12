import time
import statistics

import PARAMS
from eval import *
from utils import *

start_time = time.time()

sequence, probability_array = read_inputs(PARAMS.sequence_dir, PARAMS.confidence_dir)
index = build_index(probability_array, PARAMS.w, PARAMS.confidence_threshold)
print("Index building", time.time()-start_time)
def run():
    q, true_start_pos = gen_query_sequence(
        probability_array,
        PARAMS.query_length,
        PARAMS.substitution_rate,
        PARAMS.deletion_rate,
        PARAMS.insertion_rate,
        sanity=False,
    )
    print("True Start Position:", true_start_pos)
    candidates = find_probabilistic_matches(
        q,
        probability_array,
        index,
        PARAMS.w,
        PARAMS.score_threshold,
        PARAMS.confidence_threshold,
    )


    # Number of top x% HSPs to keep
    x_percent = 10  # For example, keep the top 10% of HSPs
    HSPs = []

    print("Starting ungapped extension")
    for candidate_start, _ in candidates:
        extension_score = ungapped_extension_probabilistic(
            q, probability_array, candidate_start, PARAMS.w, PARAMS.gap_penalty
        )
        # if extension_score >= HSP_score_threshold:
        HSPs.append((candidate_start, extension_score))

    # Sort HSPs by extension score in descending order
    HSPs = sorted(HSPs, key=lambda x: x[1], reverse=True)

    # top_x_count = int(len(HSPs) * x_percent / 100)
    top_x_count = 20
    HSPs_filtered = HSPs[:top_x_count]
    HSPs = HSPs_filtered

    print(HSPs)

    gapped_results = []
    for candidate_start, ungapped_score in HSPs:
        gapped_score, aligned_query, aligned_reference = gapped_extension(
            q, probability_array, candidate_start, PARAMS.gap_penalty
        )
        gapped_results.append(
            (gapped_score, candidate_start, aligned_query, aligned_reference)
        )
    gapped_results = sorted(gapped_results, key=lambda x: x[0], reverse=True)

    # Output results sorted by gapped score
    count = 1
    for (
        gapped_score,
        candidate_start,
        aligned_query,
        aligned_reference,
    ) in gapped_results:
        print(f"HSP{count}")
        print(f"Candidate start: {candidate_start}, Gapped score: {gapped_score}")
        print(f"Query Alignment    : {aligned_query}")
        print(f"Reference Alignment: {aligned_reference}")
        count += 1
        break

    execution_time = time.time() - start_time
    print(f"\n\nExec time: {execution_time:.4f} seconds")

    return aligned_query, aligned_reference, candidate_start, true_start_pos


correct = 0
ngrams = []
lcs = []
entropy = []
runtime = []

for _ in range(5):
    # start_time = time.time()
    top_q, top_r, candidate_start, true_start_pos = run()
    print(time.time-start_time)
    break

    if candidate_start < true_start_pos + PARAMS.query_length and candidate_start > true_start_pos - PARAMS.query_length:
        correct += 1

    score = evaluate_alignment(top_q, top_r, method="ngram", k=5)
    ngrams.append(score)

    score_lcs = evaluate_alignment(top_q, top_r, method="lcs")
    lcs.append(score_lcs)

    entropy.append(alignment_entropy(top_q, top_r))
    runtime.append(time.time()-start_time)


# print("Num correct", correct)
# mean = statistics.mean(ngrams)
# std_dev = statistics.stdev(ngrams)
# print("Ngram:", mean, "+/-", std_dev)
# mean = statistics.mean(lcs)
# std_dev = statistics.stdev(lcs)
# print("lcs:", mean, "+/-", std_dev)
# mean = statistics.mean(entropy)
# std_dev = statistics.stdev(entropy)
# print("entropy:", mean, "+/-", std_dev)
# mean = statistics.mean(runtime)
# std_dev = statistics.stdev(runtime)
# print("runtime:", mean, "+/-", std_dev)
