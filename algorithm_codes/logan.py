import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast, get_species_in_null_vector

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

pd.set_option('display.max_columns', None)

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_logan = pd.read_csv("../mechanisms/logan/log81_EdgeList_withN2.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_logan, Sp_sparse_logan, Svv_sparse_logan = create_sparse(edge_list_logan)
# Remove any 0 columns unrelated to merging
Svv_sparse_logan, Sp_sparse_logan, Sr_sparse_logan, init_col_del_logan = del_zero_col(Svv_sparse_logan, Sp_sparse_logan, Sr_sparse_logan)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_logan_np = np.array(Svv_sparse_logan, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_logan_np = np.linalg.matrix_rank(Svv_logan_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_logan = Svv_sparse_logan.shape[0] - rank_Svv_logan_np


# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_logan = create_coproduction(Sr_sparse_logan)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_logan = create_symbols(coproduction_cols_logan)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_logan, col_del_logan = merge_coprod(Sr_sparse_logan, Sp_sparse_logan, symbol_dict_logan, coproduction_cols_logan, Svv_sparse_logan)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_logan = S_merge_logan.shape[1] - Svv_sparse_logan.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants
rank_list_logan = linalg_experiment(S_merge_logan, num_experiments)

# Use SymPy built-in function to obtain left null space of S_merge_logan
N = S_merge_logan.T.nullspace()

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_logan)
rank_list_logan_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])

# Loop over the edge list to create a species index mapping dictionary 
# Species can either be in from or to columns
# If it begins with "R" then it is a reaction number, else it is a species
species_index = {}
for index, row in edge_list_logan.iterrows():
    from_node = row['from']
    to_node = row['to']
    if not from_node.startswith("R"):
        if from_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = from_node
    if not to_node.startswith("R"):
        if to_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = to_node

# Obtain the second vector in the left null space of S_merge_logan
first_null_vector = N[1]
species_in_first_null_vector = get_species_in_null_vector(first_null_vector, species_index)
print(species_in_first_null_vector)