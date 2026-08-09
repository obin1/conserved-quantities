import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("Isoprene Full")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_iso_full = pd.read_csv("../mechanisms/isoprene_full_v5/isoprene_full_v5_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_iso_full, Sp_sparse_iso_full, Svv_sparse_iso_full = create_sparse(edge_list_iso_full)
# Remove any 0 columns unrelated to merging
Svv_sparse_iso_full, Sp_sparse_iso_full, Sr_sparse_iso_full, init_col_del_iso_full = del_zero_col(Svv_sparse_iso_full, Sp_sparse_iso_full, Sr_sparse_iso_full)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_iso_full_np = np.array(Svv_sparse_iso_full, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_iso_full_np = np.linalg.matrix_rank(Svv_iso_full_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_iso_full = Svv_sparse_iso_full.shape[0] - rank_Svv_iso_full_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_iso_full = create_coproduction(Sr_sparse_iso_full)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_iso_full = create_symbols(coproduction_cols_iso_full, init_col_del_iso_full)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_iso_full, col_del_iso_full = merge_coprod(Sr_sparse_iso_full, Sp_sparse_iso_full, symbol_dict_iso_full, coproduction_cols_iso_full, Svv_sparse_iso_full)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_iso_full = S_merge_iso_full.shape[1] - Svv_sparse_iso_full.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_iso_full = linalg_experiment(S_merge_iso_full, num_experiments)
