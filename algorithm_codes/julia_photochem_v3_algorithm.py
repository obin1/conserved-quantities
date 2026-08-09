import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

print("JPM v3")

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Define stoichiometric matrix
# no Oxygen tracked in reaction 3

# new reactions: 
# reaction 14 is PAN + OH -> HCHO + CO + NO2 + H2 + O2
# reaction 15 is O3 + OH -> HO2 + O2
Svv_JPM3_hardcoded = sp.Matrix([
               [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # O₃
               [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1,  1,  0], # NO₂
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # HCHO 
               [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0,  0,  1], # HO₂  
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0, -1, -1], # OH
               [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0], # CO   
               [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # H₂   
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1,  0,  0], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  0], # PAN
               [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0], # H₂O
               [-1,  1,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0,  1,  1]])# O₂  
# Convert to sparse matrix
Svv_sparse_JPM3_hardcoded = sp.SparseMatrix(Svv_JPM3_hardcoded)

# Define product matrix
Sp_JPM3_hardcoded = sp.Matrix([
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # O₃
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # NO
               [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0], # NO₂
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], # HCHO
               [0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1], # HO₂
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # H₂O₂
               [0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0], # OH
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], # HNO₃
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0], # CO
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], # H₂
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # ALD2
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # MGLY
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0], # MCO₃
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], # PAN
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], # H₂O
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])# O₂
# Convert to sparse matrix
Sp_sparse_JPM3_hardcoded = sp.SparseMatrix(Sp_JPM3_hardcoded)

# Define reactant matrix
Sr_JPM3_hardcoded = sp.Matrix([
               [ 0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # O₃
               [ 0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0, -1,  0,  0,  0], # NO₂
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # HCHO
               [ 0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # HO₂
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,  0,  0, -1,  0, -1,  0, -1, -1,  0,  0,  0, -1, -1], # OH
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0], # PAN
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂O
               [-1,  0,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0,  0,  0]])# O₂
# Convert to sparse matrix
Sr_sparse_JPM3_hardcoded = sp.SparseMatrix(Sr_JPM3_hardcoded)

edge_list_JPM3 = pd.read_csv("../mechanisms/jpm_v3/jpm_v3_EdgeList.csv", comment="!")
# print("creating sparse matrices...")
Sr_sparse_JPM3, Sp_sparse_JPM3, Svv_sparse_JPM3 = create_sparse(edge_list_JPM3)

# Compute dimension of left null space
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_np_JPM3 = np.array(Svv_sparse_JPM3, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_JPM3_np = np.linalg.matrix_rank(Svv_np_JPM3)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_JPM3 = Svv_sparse_JPM3.shape[0] - rank_Svv_JPM3_np

# for Sr get rid of the -2 in the O2 column for reaction 3, since it does not affect the mass action law
Sr_sparse_JPM3[15, 2] = 0
Svv_sparse_JPM3[15, 2] = 0

# # assert that the hardcoded matrices are equal to the ones generated from the edge list
assert Sr_sparse_JPM3_hardcoded == Sr_sparse_JPM3, "Sr matrices are not equal"
assert Sp_sparse_JPM3_hardcoded == Sp_sparse_JPM3, "Sp matrices are not equal"
assert Svv_sparse_JPM3_hardcoded == Svv_sparse_JPM3, "Svv matrices are not equal"

init_col_del_jpm3 = []
# Identify sets of coproducing reactions
print("identifying coproduction columns...")
# Create the dictionary of symbols for coproducing groups
coproduction_cols_JPM3 = create_coproduction(Sr_sparse_JPM3)
print("creating symbolic dictionary...")
symbol_dict_JPM3 = create_symbols(coproduction_cols_JPM3, init_col_del_jpm3)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_JPM3, col_del_JPM3 = merge_coprod(Sr_sparse_JPM3, Sp_sparse_JPM3, symbol_dict_JPM3, coproduction_cols_JPM3, Svv_sparse_JPM3)
# Reinstate the O2 that was not trackedin the beginning for reaction 3
# We include the symbol for partial kinetic reaction rates after merging
S_merge_JPM3[15, 2] = -2*symbol_dict_JPM3[3]
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_JPM3 = S_merge_JPM3.shape[1] - Svv_sparse_JPM3.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge_d.shape[0] - rank(S_merge) - # stoichiometric invariants
rank_list_JPM3 = linalg_experiment(S_merge_JPM3, num_experiments)
