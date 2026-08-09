import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("Isoprene Reduced")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_iso_red = pd.read_csv("../mechanisms/isoprene_reduced_plus_v5/isoprene_reduced_plus_v5_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_iso_red, Sp_sparse_iso_red, Svv_sparse_iso_red = create_sparse(edge_list_iso_red)
# Remove any 0 columns unrelated to merging
Svv_sparse_iso_red, Sp_sparse_iso_red, Sr_sparse_iso_red, init_col_del_iso_red = del_zero_col(Svv_sparse_iso_red, Sp_sparse_iso_red, Sr_sparse_iso_red)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_iso_red_np = np.array(Svv_sparse_iso_red, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_iso_red_np = np.linalg.matrix_rank(Svv_iso_red_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_iso_red = Svv_sparse_iso_red.shape[0] - rank_Svv_iso_red_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_iso_red = create_coproduction(Sr_sparse_iso_red)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_iso_red = create_symbols(coproduction_cols_iso_red, init_col_del_iso_red)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_iso_red, col_del_iso_red = merge_coprod(Sr_sparse_iso_red, Sp_sparse_iso_red, symbol_dict_iso_red, coproduction_cols_iso_red, Svv_sparse_iso_red)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_iso_red = S_merge_iso_red.shape[1] - Svv_sparse_iso_red.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_iso_red = linalg_experiment(S_merge_iso_red, num_experiments)
