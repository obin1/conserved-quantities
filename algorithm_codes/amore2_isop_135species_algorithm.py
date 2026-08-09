import pandas as pd
import numpy as np
from sympy import S
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

print("Amore")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_amore = pd.read_csv("../mechanisms/amore2_isop_135species/amore2_isop_135species_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_amore, Sp_sparse_amore, Svv_sparse_amore = create_sparse(edge_list_amore)
# Remove any 0 columns unrelated to merging
Svv_sparse_amore, Sp_sparse_amore, Sr_sparse_amore, init_col_del_amore = del_zero_col(Svv_sparse_amore, Sp_sparse_amore, Sr_sparse_amore)

# Remove rows of species that do not participate in any reactions
zero_rows_indices = []
for i in range(Svv_sparse_amore.rows):
    row = Svv_sparse_amore.row(i)
    if all(element == S.Zero for element in row):
        zero_rows_indices.append(i)
for i in reversed(zero_rows_indices):
    Svv_sparse_amore.row_del(i)
    Sr_sparse_amore.row_del(i)
    Sp_sparse_amore.row_del(i)

print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_amore_np = np.array(Svv_sparse_amore, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_amore_np = np.linalg.matrix_rank(Svv_amore_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_amore = Svv_sparse_amore.shape[0] - rank_Svv_amore_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_amore = create_coproduction(Sr_sparse_amore)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_amore = create_symbols(coproduction_cols_amore, init_col_del_amore)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_amore, col_del_amore = merge_coprod(Sr_sparse_amore, Sp_sparse_amore, symbol_dict_amore, coproduction_cols_amore, Svv_sparse_amore)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_amore = S_merge_amore.shape[1] - Svv_sparse_amore.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_amore = linalg_experiment(S_merge_amore, num_experiments)
