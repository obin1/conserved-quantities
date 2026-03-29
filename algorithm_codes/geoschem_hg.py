import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_gchg= pd.read_csv("../mechanisms/geoschem-hg/Hg_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_gchg, Sp_sparse_gchg, Svv_sparse_gchg= create_sparse(edge_list_gchg)
# Remove any 0 columns unrelated to merging
Svv_sparse_gchg, init_col_del_gchg = del_zero_col(Svv_sparse_gchg)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_gchg_np = np.array(Svv_sparse_gchg, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_gchg_np = np.linalg.matrix_rank(Svv_gchg_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_gchg= Svv_sparse_gchg.shape[0] - rank_Svv_gchg_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_gchg= create_coproduction(Sr_sparse_gchg)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_gchg= create_symbols(coproduction_cols_gchg)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_gchg, col_del_gchg= merge_coprod(Sr_sparse_gchg, Sp_sparse_gchg, symbol_dict_gchg, coproduction_cols_gchg, Svv_sparse_gchg)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_gchg= S_merge_gchg.shape[1] - Svv_sparse_gchg.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_gchg= linalg_experiment(S_merge_gchg, num_experiments)

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_gchg)
rank_list_gchg_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])

