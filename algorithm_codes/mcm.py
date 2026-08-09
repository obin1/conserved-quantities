import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, linalg_experiment_fast, numeric_sparse_matrix_fast_combined

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 1

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_mcm = pd.read_csv("../mechanisms/mcm_v3.3.1/mcm_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_mcm, Sp_sparse_mcm, Svv_sparse_mcm = create_sparse(edge_list_mcm)
# Remove any 0 columns unrelated to merging
Svv_sparse_mcm, Sp_sparse_mcm, Sr_sparse_mcm, init_col_del_mcm = del_zero_col(Svv_sparse_mcm, Sp_sparse_mcm, Sr_sparse_mcm)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_mcm_np = np.array(Svv_sparse_mcm, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_mcm_np = np.linalg.matrix_rank(Svv_mcm_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_mcm = Svv_sparse_mcm.shape[0] - rank_Svv_mcm_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_mcm = create_coproduction(Sr_sparse_mcm)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_mcm = create_symbols(coproduction_cols_mcm, init_col_del_mcm)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_mcm, col_del_mcm = merge_coprod(Sr_sparse_mcm, Sp_sparse_mcm, symbol_dict_mcm, coproduction_cols_mcm, Svv_sparse_mcm)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_mcm = S_merge_mcm.shape[1] - Svv_sparse_mcm.shape[1]

# Convert S_merge (SymPy sparse matrix) into SciPy Compressed Sparse Column format matrix
A = numeric_sparse_matrix_fast_combined(S_merge_mcm)
# Project S_merge to a smaller dimension then perform linear algebra experimentation
rank_list_mcm = linalg_experiment_fast(A, num_experiments, A.shape[0])
