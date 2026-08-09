import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("JAM")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_jam = pd.read_csv("../mechanisms/jam/jam_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_jam, Sp_sparse_jam, Svv_sparse_jam = create_sparse(edge_list_jam)
# Remove any 0 columns unrelated to merging
Svv_sparse_jam, Sp_sparse_jam, Sr_sparse_jam, init_col_del_jam = del_zero_col(Svv_sparse_jam, Sp_sparse_jam, Sr_sparse_jam)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_jam_np = np.array(Svv_sparse_jam, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_jam_np = np.linalg.matrix_rank(Svv_jam_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_jam = Svv_sparse_jam.shape[0] - rank_Svv_jam_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_jam = create_coproduction(Sr_sparse_jam)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_jam = create_symbols(coproduction_cols_jam, init_col_del_jam)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_jam, col_del_jam = merge_coprod(Sr_sparse_jam, Sp_sparse_jam, symbol_dict_jam, coproduction_cols_jam, Svv_sparse_jam)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_jam = S_merge_jam.shape[1] - Svv_sparse_jam.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_jam = linalg_experiment(S_merge_jam, num_experiments)

