import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_gc13_3_0 = pd.read_csv("../mechanisms/geoschem-13.3.0/gckpp_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_gc13_3_0, Sp_sparse_gc13_3_0, Svv_sparse_gc13_3_0 = create_sparse(edge_list_gc13_3_0)
# Remove any 0 columns unrelated to merging
Svv_sparse_gc13_3_0, Sp_sparse_gc13_3_0, Sr_sparse_gc13_3_0, init_col_del_gc13_3_0 = del_zero_col(Svv_sparse_gc13_3_0, Sp_sparse_gc13_3_0, Sr_sparse_gc13_3_0)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_gc13_3_0_np = np.array(Svv_sparse_gc13_3_0, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_gc13_3_0_np = np.linalg.matrix_rank(Svv_gc13_3_0_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_gc13_3_0 = Svv_sparse_gc13_3_0.shape[0] - rank_Svv_gc13_3_0_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_gc13_3_0 = create_coproduction(Sr_sparse_gc13_3_0)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_gc13_3_0 = create_symbols(coproduction_cols_gc13_3_0)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_gc13_3_0, col_del_gc13_3_0 = merge_coprod(Sr_sparse_gc13_3_0, Sp_sparse_gc13_3_0, symbol_dict_gc13_3_0, coproduction_cols_gc13_3_0, Svv_sparse_gc13_3_0)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_gc13_3_0 = S_merge_gc13_3_0.shape[1] - Svv_sparse_gc13_3_0.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_gc13_3_0 = linalg_experiment(S_merge_gc13_3_0, num_experiments)

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_gc13_3_0)
rank_list_gc13_3_0_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])