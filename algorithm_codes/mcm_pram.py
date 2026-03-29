import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, linalg_experiment_fast, numeric_sparse_matrix_fast_combined

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 1

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_pram = pd.read_csv("../mechanisms/mcm-pram/mcm_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_pram, Sp_sparse_pram, Svv_sparse_pram = create_sparse(edge_list_pram)
# Remove any 0 columns unrelated to merging
Svv_sparse_pram, init_col_del_pram = del_zero_col(Svv_sparse_pram)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_pram_np = np.array(Svv_sparse_pram, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_pram_np = np.linalg.matrix_rank(Svv_pram_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_pram = Svv_sparse_pram.shape[0] - rank_Svv_pram_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_pram = create_coproduction(Sr_sparse_pram)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_pram = create_symbols(coproduction_cols_pram)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_pram, col_del_pram = merge_coprod(Sr_sparse_pram, Sp_sparse_pram, symbol_dict_pram, coproduction_cols_pram, Svv_sparse_pram)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_pram = S_merge_pram.shape[1] - Svv_sparse_pram.shape[1]

# Convert S_merge (SymPy sparse matrix) into SciPy Compressed Sparse Column format matrix
A = numeric_sparse_matrix_fast_combined(S_merge_pram)
# Project S_merge to a smaller dimension then perform linear algebra experimentation
rank_list_pram = linalg_experiment_fast(A, num_experiments, A.shape[0])
