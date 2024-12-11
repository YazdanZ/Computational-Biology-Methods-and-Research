# ----- Inputs ------ #
sequence_dir = "chr22.maf.ancestors.42000000.complete.boreo.fa"
confidence_dir = "chr22.maf.ancestors.42000000.complete.boreo.conf"
rows_to_nucleotides = {0: "A", 1: "C", 2: "G", 3: "T"}
nucleotides_to_rows = {"A": 0, "C": 1, "G": 2, "T": 3}

# ----- Indexing ------ #
w = 12  # the w-mer
threshold = 0.90  # threshold used to create w-mer index
# HSP_score_threshold = -4


# ----- Query Generation ------ #
query_length = 100
substitution_rate = 2  # in %
deletion_rate = 2  # in %
insertion_rate = 2  # in %

# ----- scoring ------ #
score_threshold = 0.5
confidence_threshold = 0.6
gap_penalty = 0.5
