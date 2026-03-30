import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_gc = pd.read_csv("../mechanisms/geos-chem-v14/gckpp_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_gc, Sp_sparse_gc, Svv_sparse_gc = create_sparse(edge_list_gc)
# Remove any 0 columns unrelated to merging
Svv_sparse_gc, Sp_sparse_gc, Sr_sparse_gc, init_col_del_gc = del_zero_col(Svv_sparse_gc, Sp_sparse_gc, Sr_sparse_gc)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_gc_np = np.array(Svv_sparse_gc, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_gc_np = np.linalg.matrix_rank(Svv_gc_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_gc = Svv_gc_np.shape[0] - rank_Svv_gc_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_gc = create_coproduction(Sr_sparse_gc)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_gc = create_symbols(coproduction_cols_gc)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_gc, col_del_gc = merge_coprod(Sr_sparse_gc, Sp_sparse_gc, symbol_dict_gc, coproduction_cols_gc, Svv_sparse_gc)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_gc = S_merge_gc.shape[1] - Svv_sparse_gc.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_gc = linalg_experiment(S_merge_gc, num_experiments)
