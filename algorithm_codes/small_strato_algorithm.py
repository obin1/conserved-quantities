import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_small_strato = pd.read_csv("../mechanisms/small_strato/small_strato_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_small_strato, Sp_sparse_small_strato, Svv_sparse_small_strato = create_sparse(edge_list_small_strato)
# Remove any 0 columns unrelated to merging
Svv_sparse_small_strato, Sp_sparse_small_strato, Sr_sparse_small_strato, init_col_del_small_strato = del_zero_col(Svv_sparse_small_strato, Sp_sparse_small_strato, Sr_sparse_small_strato)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_small_strato_np = np.array(Svv_sparse_small_strato, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_small_strato_np = np.linalg.matrix_rank(Svv_small_strato_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_small_strato = Svv_sparse_small_strato.shape[0] - rank_Svv_small_strato_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_small_strato = create_coproduction(Sr_sparse_small_strato)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_small_strato = create_symbols(coproduction_cols_small_strato)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_small_strato, col_del_small_strato = merge_coprod(Sr_sparse_small_strato, Sp_sparse_small_strato, symbol_dict_small_strato, coproduction_cols_small_strato, Svv_sparse_small_strato)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_small_strato = S_merge_small_strato.shape[1] - Svv_sparse_small_strato.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_small_strato = linalg_experiment(S_merge_small_strato, num_experiments)