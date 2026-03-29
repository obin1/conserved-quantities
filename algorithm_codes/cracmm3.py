# WITH DELTA_C, DELTA_N, DELTA_SI
import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cracmm3 = pd.read_csv("../mechanisms/cracmm3/cracmm3_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cracmm3, Sp_sparse_cracmm3, Svv_sparse_cracmm3 = create_sparse(edge_list_cracmm3)
# Remove any 0 columns unrelated to merging
Svv_sparse_cracmm3, init_col_del_cracmm3 = del_zero_col(Svv_sparse_cracmm3)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cracmm3_np = np.array(Svv_sparse_cracmm3, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cracmm3_np = np.linalg.matrix_rank(Svv_cracmm3_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cracmm3 = Svv_sparse_cracmm3.shape[0] - rank_Svv_cracmm3_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cracmm3 = create_coproduction(Sr_sparse_cracmm3)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cracmm3 = create_symbols(coproduction_cols_cracmm3)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cracmm3, col_del_cracmm3 = merge_coprod(Sr_sparse_cracmm3, Sp_sparse_cracmm3, symbol_dict_cracmm3, coproduction_cols_cracmm3, Svv_sparse_cracmm3)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cracmm3 = S_merge_cracmm3.shape[1] - Svv_sparse_cracmm3.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_cracmm3 = linalg_experiment(S_merge_cracmm3, num_experiments)


#%%
# WITHOUT DELTA_C, DELTA_N, DELTA_SI

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_cracmm3e = pd.read_csv("../mechanisms/cracmm3/cracmm3_EdgeList_eliminated.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_cracmm3e, Sp_sparse_cracmm3e, Svv_sparse_cracmm3e = create_sparse(edge_list_cracmm3e)
# Remove any 0 columns unrelated to merging
Svv_sparse_cracmm3e, init_col_del_cracmm3e = del_zero_col(Svv_sparse_cracmm3e)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_cracmm3e_np = np.array(Svv_sparse_cracmm3e, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_cracmm3e_np = np.linalg.matrix_rank(Svv_cracmm3e_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_cracmm3e = Svv_sparse_cracmm3e.shape[0] - rank_Svv_cracmm3e_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_cracmm3e = create_coproduction(Sr_sparse_cracmm3e)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_cracmm3e = create_symbols(coproduction_cols_cracmm3e)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_cracmm3e, col_del_cracmm3e = merge_coprod(Sr_sparse_cracmm3e, Sp_sparse_cracmm3e, symbol_dict_cracmm3e, coproduction_cols_cracmm3e, Svv_sparse_cracmm3e)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_cracmm3e = S_merge_cracmm3e.shape[1] - Svv_sparse_cracmm3e.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants    
rank_list_cracmm3e = linalg_experiment(S_merge_cracmm3e, num_experiments)
