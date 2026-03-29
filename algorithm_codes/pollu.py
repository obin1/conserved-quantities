import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_pollu = pd.read_csv("../mechanisms/pollu/pollu_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_pollu, Sp_sparse_pollu, Svv_sparse_pollu = create_sparse(edge_list_pollu)
# Remove any 0 columns unrelated to merging
Svv_sparse_pollu, init_col_del_pollu = del_zero_col(Svv_sparse_pollu)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_pollu_np = np.array(Svv_sparse_pollu, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_pollu_np = np.linalg.matrix_rank(Svv_pollu_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_pollu = Svv_sparse_pollu.shape[0] - rank_Svv_pollu_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_pollu = create_coproduction(Sr_sparse_pollu)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_pollu = create_symbols(coproduction_cols_pollu)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_pollu, col_del_pollu = merge_coprod(Sr_sparse_pollu, Sp_sparse_pollu, symbol_dict_pollu, coproduction_cols_pollu, Svv_sparse_pollu)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_pollu = S_merge_pollu.shape[1] - Svv_sparse_pollu.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_pollu = linalg_experiment(S_merge_pollu, num_experiments)
    