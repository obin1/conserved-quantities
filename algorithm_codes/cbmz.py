import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cbmz = pd.read_csv("../mechanisms/cbmz/CBMZ_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cbmz, Sp_sparse_cbmz, Svv_sparse_cbmz = create_sparse(edge_list_cbmz)
# Remove any 0 columns unrelated to merging
Svv_sparse_cbmz, Sp_sparse_cbmz, Sr_sparse_cbmz, init_col_del_cbmz = del_zero_col(Svv_sparse_cbmz, Sp_sparse_cbmz, Sr_sparse_cbmz)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cbmz_np = np.array(Svv_sparse_cbmz, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cbmz_np = np.linalg.matrix_rank(Svv_cbmz_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cbmz = Svv_sparse_cbmz.shape[0] - rank_Svv_cbmz_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cbmz = create_coproduction(Sr_sparse_cbmz)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cbmz = create_symbols(coproduction_cols_cbmz)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cbmz, col_del_cbmz = merge_coprod(Sr_sparse_cbmz, Sp_sparse_cbmz, symbol_dict_cbmz, coproduction_cols_cbmz, Svv_sparse_cbmz)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cbmz = S_merge_cbmz.shape[1] - Svv_sparse_cbmz.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_cbmz = linalg_experiment(S_merge_cbmz, num_experiments)