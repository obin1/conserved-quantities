import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_radm2 = pd.read_csv("../mechanisms/radm2/RADM2_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_radm2, Sp_sparse_radm2, Svv_sparse_radm2 = create_sparse(edge_list_radm2)
# Remove any 0 columns unrelated to merging
Svv_sparse_radm2, Sp_sparse_radm2, Sr_sparse_radm2, init_col_del_radm2 = del_zero_col(Svv_sparse_radm2, Sp_sparse_radm2, Sr_sparse_radm2)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_radm2_np = np.array(Svv_sparse_radm2, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_radm2_np = np.linalg.matrix_rank(Svv_radm2_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_radm2 = Svv_sparse_radm2.shape[0] - rank_Svv_radm2_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_radm2 = create_coproduction(Sr_sparse_radm2)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_radm2 = create_symbols(coproduction_cols_radm2)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_radm2, col_del_radm2 = merge_coprod(Sr_sparse_radm2, Sp_sparse_radm2, symbol_dict_radm2, coproduction_cols_radm2, Svv_sparse_radm2)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_radm2 = S_merge_radm2.shape[1] - Svv_sparse_radm2.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_radm2 = linalg_experiment(S_merge_radm2, num_experiments)
