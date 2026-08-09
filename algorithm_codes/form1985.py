import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("Form1985")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_form1985 = pd.read_csv("../mechanisms/form1985/form1985_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_form1985, Sp_sparse_form1985, Svv_sparse_form1985 = create_sparse(edge_list_form1985)
# Remove any 0 columns unrelated to merging
Svv_sparse_form1985, Sp_sparse_form1985, Sr_sparse_form1985, init_col_del_form1985 = del_zero_col(Svv_sparse_form1985, Sp_sparse_form1985, Sr_sparse_form1985)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_form1985_np = np.array(Svv_sparse_form1985, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_form1985_np = np.linalg.matrix_rank(Svv_form1985_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_form1985 = Svv_sparse_form1985.shape[0] - rank_Svv_form1985_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_form1985 = create_coproduction(Sr_sparse_form1985)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_form1985 = create_symbols(coproduction_cols_form1985, init_col_del_form1985)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_form1985, col_del_form1985 = merge_coprod(Sr_sparse_form1985, Sp_sparse_form1985, symbol_dict_form1985, coproduction_cols_form1985, Svv_sparse_form1985)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_form1985 = S_merge_form1985.shape[1] - Svv_sparse_form1985.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_form1985 = linalg_experiment(S_merge_form1985, num_experiments)

