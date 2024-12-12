import math 

def alignment_entropy(query, reference):
  """Calculates the entropy score of an alignment.

  Args:
    query: The query sequence.
    reference: The reference sequence.

  Returns:
    The entropy score of the alignment.
  """

  if len(query) != len(reference):
    raise ValueError("Query and reference sequences must have the same length.")

  # Create a dictionary to store the frequency of character pairs
  pair_counts = {}

  for q_char, r_char in zip(query, reference):
    pair = (q_char, r_char)
    pair_counts[pair] = pair_counts.get(pair, 0) + 1

  # Calculate the total number of pairs
  total_pairs = sum(pair_counts.values())

  # Calculate the entropy
  entropy = 0
  for count in pair_counts.values():
    prob = count / total_pairs
    entropy -= prob * math.log2(prob)

  return entropy

def ngram_similarity(query, reference, k=5):
    """
    Computes the n-gram overlap similarity between two sequences (query and reference).
    
    :param query: The query sequence string.
    :param reference: The reference sequence string.
    :param k: The length of the n-grams to consider. Default is 5.
    
    :return: A similarity score (Jaccard Index) between the two sequences.
    """
    # Generate k-mers (substrings of length k)
    def generate_kmers(sequence, k):
        return set(sequence[i:i+k] for i in range(len(sequence) - k + 1))
    
    query_kmers = generate_kmers(query, k)
    reference_kmers = generate_kmers(reference, k)
    
    # Calculate the intersection and union of the two sets of k-mers
    intersection = query_kmers.intersection(reference_kmers)
    union = query_kmers.union(reference_kmers)
    
    # Compute the Jaccard index (intersection over union)
    if len(union) == 0:
        return 0.0  # Avoid division by zero if no common k-mers exist
    
    similarity_score = len(intersection) / len(union)
    return similarity_score


def lcs_similarity(query, reference):
    """
    Computes the Longest Common Subsequence (LCS) similarity score.
    
    :param query: The query sequence string.
    :param reference: The reference sequence string.
    
    :return: A similarity score based on the length of the LCS.
    """
    # Initialize the DP table
    m, n = len(query), len(reference)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if query[i - 1] == reference[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    
    lcs_length = dp[m][n]
    
    max_length = max(m, n)
    lcs_score = lcs_length / max_length
    return lcs_score


def evaluate_alignment(query, reference, method="ngram", k=5):
    if method == "ngram":
        return ngram_similarity(query, reference, k)
    elif method == "lcs":
        return lcs_similarity(query, reference)
    else:
        raise ValueError("Invalid method. Choose either 'ngram' or 'lcs'.")

