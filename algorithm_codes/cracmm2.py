import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cracmm = pd.read_csv("../mechanisms/cracmm2/cracmm2_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cracmm, Sp_sparse_cracmm, Svv_sparse_cracmm = create_sparse(edge_list_cracmm)
# Remove any 0 columns unrelated to merging
Svv_sparse_cracmm, Sp_sparse_cracmm, Sr_sparse_cracmm, init_col_del_cracmm = del_zero_col(Svv_sparse_cracmm, Sp_sparse_cracmm, Sr_sparse_cracmm)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cracmm_np = np.array(Svv_sparse_cracmm, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cracmm_np = np.linalg.matrix_rank(Svv_cracmm_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cracmm = Svv_sparse_cracmm.shape[0] - rank_Svv_cracmm_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cracmm = create_coproduction(Sr_sparse_cracmm)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cracmm = create_symbols(coproduction_cols_cracmm)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cracmm, col_del_cracmm = merge_coprod(Sr_sparse_cracmm, Sp_sparse_cracmm, symbol_dict_cracmm, coproduction_cols_cracmm, Svv_sparse_cracmm)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cracmm = S_merge_cracmm.shape[1] - Svv_sparse_cracmm.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_cracmm = linalg_experiment(S_merge_cracmm, num_experiments)
