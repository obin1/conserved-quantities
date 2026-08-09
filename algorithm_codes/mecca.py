import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("MECCA")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 1

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_mecca = pd.read_csv("../mechanisms/mecca_v4.6.0/gas_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_mecca, Sp_sparse_mecca, Svv_sparse_mecca = create_sparse(edge_list_mecca)
# Remove any 0 columns unrelated to merging
Svv_sparse_mecca, Sp_sparse_mecca, Sr_sparse_mecca, init_col_del_mecca = del_zero_col(Svv_sparse_mecca, Sp_sparse_mecca, Sr_sparse_mecca)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_mecca_np = np.array(Svv_sparse_mecca, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_mecca_np = np.linalg.matrix_rank(Svv_mecca_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_mecca = Svv_sparse_mecca.shape[0] - rank_Svv_mecca_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_mecca = create_coproduction(Sr_sparse_mecca)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_mecca = create_symbols(coproduction_cols_mecca, init_col_del_mecca)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_mecca, col_del_mecca = merge_coprod(Sr_sparse_mecca, Sp_sparse_mecca, symbol_dict_mecca, coproduction_cols_mecca, Svv_sparse_mecca)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_mecca = S_merge_mecca.shape[1] - Svv_sparse_mecca.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_mecca = linalg_experiment(S_merge_mecca, num_experiments)

