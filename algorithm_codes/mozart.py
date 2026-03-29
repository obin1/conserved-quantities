import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_mozart = pd.read_csv("../mechanisms/mozart/MOZART_4_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_mozart, Sp_sparse_mozart, Svv_sparse_mozart = create_sparse(edge_list_mozart)
# Remove any 0 columns unrelated to merging
Svv_sparse_mozart, init_col_del_mozart = del_zero_col(Svv_sparse_mozart)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_mozart_np = np.array(Svv_sparse_mozart, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_mozart_np = np.linalg.matrix_rank(Svv_mozart_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_mozart = Svv_sparse_mozart.shape[0] - rank_Svv_mozart_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_mozart = create_coproduction(Sr_sparse_mozart)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_mozart = create_symbols(coproduction_cols_mozart)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_mozart, col_del_mozart = merge_coprod(Sr_sparse_mozart, Sp_sparse_mozart, symbol_dict_mozart, coproduction_cols_mozart, Svv_sparse_mozart)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_mozart = S_merge_mozart.shape[1] - Svv_sparse_mozart.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_mozart = linalg_experiment(S_merge_mozart, num_experiments)
  