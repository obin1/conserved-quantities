import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Define stoichiometric matrix
# no Oxygen tracked in reaction 3
Svv_JPM2 = sp.Matrix([
               [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O₃
               [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1], # NO₂
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0], # HCHO 
               [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0], # HO₂  
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0], # OH
               [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0], # CO   
               [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1], # PAN
               [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0], # H₂O
               [-1,  1,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0]])# O₂  
# Convert to sparse matrix
Svv_sparse_JPM2 = sp.SparseMatrix(Svv_JPM2)

# Define reactant matrix
Sp_JPM2 = sp.Matrix([
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
# Convert to sparse matrix
Sp_sparse_JPM2 = sp.SparseMatrix(Sp_JPM2)

# Define product matrix
Sr_JPM2 = sp.Matrix([
               [ 0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
               [-1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0, -1,  0],
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0],
               [ 0,  0,  0,  0, -1,  0, -1,  0, -1, -1,  0,  0,  0], 
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [-1,  0,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0]])
# Convert to sparse matrix
Sr_sparse_JPM2 = sp.SparseMatrix(Sr_JPM2)

# edge_list_JPM2 = pd.read_csv("../mechanisms/jpm_v2/jpm_v2_EdgeList.csv", comment="!")

# print("creating sparse matrices...")
# Sr_sparse_JPM2, Sp_sparse_JPM2, Svv_sparse_JPM2 = create_sparse(edge_list_JPM2)
# Svv_sparse_JPM2, init_col_del_JPM2 = del_zero_col(Svv_sparse_JPM2)

# Compute dimension of left null space in two ways 
# (SymPy nullspace operations sometimes result in a different number due to precision errors)
print("computing dimension of nullspace...")
# Compute dimension of left null space using SymPy nullspace operation
stoichiometric_invariants_JPM2 = len(Svv_sparse_JPM2.T.nullspace())
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_np_JPM2 = np.array(Svv_sparse_JPM2, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_JPM2_np = np.linalg.matrix_rank(Svv_np_JPM2)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_JPM2 = Svv_sparse_JPM2.shape[0] - rank_Svv_JPM2_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
# Create the dictionary of symbols for coproducing groups
coproduction_cols_JPM2 = create_coproduction(Sr_sparse_JPM2)
print("creating symbolic dictionary...")
symbol_dict_JPM2 = create_symbols(coproduction_cols_JPM2)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_JPM2, col_del_JPM2 = merge_coprod(Sr_sparse_JPM2, Sp_sparse_JPM2, symbol_dict_JPM2, coproduction_cols_JPM2, Svv_sparse_JPM2)
# Reinstate the O2 that was not trackedin the beginning for reaction 3
# We include the symbol for partial kinetic reaction rates after merging
S_merge_JPM2[15, 2] = -2*symbol_dict_JPM2[3]
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_jpm2 = S_merge_JPM2.shape[1] - Svv_sparse_JPM2.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge_d.shape[0] - rank(S_merge) - # stoichiometric invariants
rank_list_jpm2 = linalg_experiment(S_merge_JPM2, num_experiments)
