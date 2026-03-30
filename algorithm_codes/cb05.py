import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cb05 = pd.read_csv("../mechanisms/cb05/CB05TUCl_EPA_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cb05, Sp_sparse_cb05, Svv_sparse_cb05 = create_sparse(edge_list_cb05)
# Remove any 0 columns unrelated to merging
Svv_sparse_cb05, Sp_sparse_cb05, Sr_sparse_cb05, init_col_del_cb05 = del_zero_col(Svv_sparse_cb05, Sp_sparse_cb05, Sr_sparse_cb05)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cb05_np = np.array(Svv_sparse_cb05, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cb05_np = np.linalg.matrix_rank(Svv_cb05_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cb05 = Svv_sparse_cb05.shape[0] - rank_Svv_cb05_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cb05 = create_coproduction(Sr_sparse_cb05)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cb05 = create_symbols(coproduction_cols_cb05)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cb05, col_del_cb05 = merge_coprod(Sr_sparse_cb05, Sp_sparse_cb05, symbol_dict_cb05, coproduction_cols_cb05, Svv_sparse_cb05)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cb05 = S_merge_cb05.shape[1] - Svv_sparse_cb05.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_cb05 = linalg_experiment(S_merge_cb05, num_experiments)

# For code testing purposes:
# Use MCM functions numeric_sparse_matrix_fast_combined and linalg_experiment_fast to perform linear algebra experimentation
A = numeric_sparse_matrix_fast_combined(S_merge_cb05)
rank_list_cb05_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])
