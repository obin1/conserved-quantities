import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_pact1d = pd.read_csv("../mechanisms/pact1d/mech_EdgeList.csv", comment="!")
# this mechanism is from pact-1d halogens, 
# paper https://doi.org/10.1029/2021JD036140
# zenodo https://doi.org/10.5281/zenodo.6045999

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_pact1d, Sp_sparse_pact1d, Svv_sparse_pact1d = create_sparse(edge_list_pact1d)
# Remove any 0 columns unrelated to merging
Svv_sparse_pact1d, init_col_del_pact1d = del_zero_col(Svv_sparse_pact1d)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_pact1d_np = np.array(Svv_sparse_pact1d, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_pact1d_np = np.linalg.matrix_rank(Svv_pact1d_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_pact1d = Svv_sparse_pact1d.shape[0] - rank_Svv_pact1d_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_pact1d = create_coproduction(Sr_sparse_pact1d)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_pact1d = create_symbols(coproduction_cols_pact1d)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_pact1d, col_del_pact1d = merge_coprod(Sr_sparse_pact1d, Sp_sparse_pact1d, symbol_dict_pact1d, coproduction_cols_pact1d, Svv_sparse_pact1d)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_pact1d = S_merge_pact1d.shape[1] - Svv_sparse_pact1d.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_pact1d = linalg_experiment(S_merge_pact1d, num_experiments)
