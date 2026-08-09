import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("E3SM")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_chemuci = pd.read_csv("../mechanisms/e3sm-chem/chemuci_EdgeList.csv", comment="!")
# this mechanism is from pact-1d halogens, 
# paper https://doi.org/10.1029/2021JD036140
# zenodo https://doi.org/10.5281/zenodo.6045999

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_chemuci, Sp_sparse_chemuci, Svv_sparse_chemuci = create_sparse(edge_list_chemuci)
# Remove any 0 columns unrelated to merging
Svv_sparse_chemuci, Sp_sparse_chemuci, Sr_sparse_chemuci, init_col_del_uci = del_zero_col(Svv_sparse_chemuci, Sp_sparse_chemuci, Sr_sparse_chemuci)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_chemuci_np = np.array(Svv_sparse_chemuci, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_chemuci_np = np.linalg.matrix_rank(Svv_chemuci_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_chemuci = Svv_sparse_chemuci.shape[0] - rank_Svv_chemuci_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_chemuci = create_coproduction(Sr_sparse_chemuci)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_chemuci = create_symbols(coproduction_cols_chemuci, init_col_del_uci)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_chemuci, col_del_chemuci = merge_coprod(Sr_sparse_chemuci, Sp_sparse_chemuci, symbol_dict_chemuci, coproduction_cols_chemuci, Svv_sparse_chemuci)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_chemuci = S_merge_chemuci.shape[1] - Svv_sparse_chemuci.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_chemuci = linalg_experiment(S_merge_chemuci, num_experiments)
