import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

print("Isoprene Full - no fixed spc")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_iso_vo = pd.read_csv("../mechanisms/isoprene_full_v5_nofixedspc/isoprene_full_v5_nofixedspc_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_iso_vo, Sp_sparse_iso_vo, Svv_sparse_iso_vo = create_sparse(edge_list_iso_vo)
# Remove any 0 columns unrelated to merging
Svv_sparse_iso_vo, Sp_sparse_iso_vo, Sr_sparse_iso_vo, init_col_del_iso_vo = del_zero_col(Svv_sparse_iso_vo, Sp_sparse_iso_vo, Sr_sparse_iso_vo)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_iso_vo_np = np.array(Svv_sparse_iso_vo, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_iso_vo_np = np.linalg.matrix_rank(Svv_iso_vo_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_iso_vo = Svv_sparse_iso_vo.shape[0] - rank_Svv_iso_vo_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_iso_vo = create_coproduction(Sr_sparse_iso_vo)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_iso_vo = create_symbols(coproduction_cols_iso_vo, init_col_del_iso_vo)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_iso_vo, col_del_iso_vo = merge_coprod(Sr_sparse_iso_vo, Sp_sparse_iso_vo, symbol_dict_iso_vo, coproduction_cols_iso_vo, Svv_sparse_iso_vo)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_iso_vo = S_merge_iso_vo.shape[1] - Svv_sparse_iso_vo.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_iso_vo = linalg_experiment(S_merge_iso_vo, num_experiments)

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_iso_vo)
rank_list_iso_vo_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])


