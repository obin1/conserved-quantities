import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_saprc99 = pd.read_csv("../mechanisms/saprc99/saprc99_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_saprc99, Sp_sparse_saprc99, Svv_sparse_saprc99 = create_sparse(edge_list_saprc99)
# Remove any 0 columns unrelated to merging
Svv_sparse_saprc99, Sp_sparse_saprc99, Sr_sparse_saprc99, init_col_del_saprc99 = del_zero_col(Svv_sparse_saprc99, Sp_sparse_saprc99, Sr_sparse_saprc99)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_saprc99_np = np.array(Svv_sparse_saprc99, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_saprc99_np = np.linalg.matrix_rank(Svv_saprc99_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_saprc99 = Svv_sparse_saprc99.shape[0] - rank_Svv_saprc99_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_saprc99 = create_coproduction(Sr_sparse_saprc99)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_saprc99 = create_symbols(coproduction_cols_saprc99)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_saprc99, col_del_saprc99 = merge_coprod(Sr_sparse_saprc99, Sp_sparse_saprc99, symbol_dict_saprc99, coproduction_cols_saprc99, Svv_sparse_saprc99)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_saprc99 = S_merge_saprc99.shape[1] - Svv_sparse_saprc99.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_saprc99 = linalg_experiment(S_merge_saprc99, num_experiments)
