import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

print("CRI v2.2")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cri = pd.read_csv("../mechanisms/cri-v2.2/cri_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cri, Sp_sparse_cri, Svv_sparse_cri = create_sparse(edge_list_cri)
# Remove any 0 columns unrelated to merging
Svv_sparse_cri, Sp_sparse_cri, Sr_sparse_cri, init_col_del_cri = del_zero_col(Svv_sparse_cri, Sp_sparse_cri, Sr_sparse_cri)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cri_np = np.array(Svv_sparse_cri, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cri_np = np.linalg.matrix_rank(Svv_cri_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cri = Svv_sparse_cri.shape[0] - rank_Svv_cri_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cri = create_coproduction(Sr_sparse_cri)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cri = create_symbols(coproduction_cols_cri, init_col_del_cri)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cri, col_del_cri = merge_coprod(Sr_sparse_cri, Sp_sparse_cri, symbol_dict_cri, coproduction_cols_cri, Svv_sparse_cri)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cri = S_merge_cri.shape[1] - Svv_sparse_cri.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_cri = linalg_experiment(S_merge_cri, num_experiments)

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_cri)
rank_list_cri_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])
