import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_JPM = pd.read_csv("../mechanisms/jpm_v1/jpm_v1_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_JPM, Sp_sparse_JPM, Svv_sparse_JPM = create_sparse(edge_list_JPM)
# Remove any 0 columns unrelated to merging
Svv_sparse_JPM, init_col_del_JPM = del_zero_col(Svv_sparse_JPM)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_np_JPM = np.array(Svv_sparse_JPM, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_JPM_np = np.linalg.matrix_rank(Svv_np_JPM)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_JPM = Svv_sparse_JPM.shape[0] - rank_Svv_JPM_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_JPM = create_coproduction(Sr_sparse_JPM)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_JPM = create_symbols(coproduction_cols_JPM)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_JPM, col_del_JPM = merge_coprod(Sr_sparse_JPM, Sp_sparse_JPM, symbol_dict_JPM, coproduction_cols_JPM, Svv_sparse_JPM)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_jpm = S_merge_JPM.shape[1] - Svv_sparse_JPM.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_jpm = linalg_experiment(S_merge_JPM, num_experiments)
