import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

SEED = 42
np.random.seed(SEED)

num_experiments = 10

# # no Oxygen tracked in reaction 3
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
Svv_sparse_JPM2 = sp.SparseMatrix(Svv_JPM2)

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
Sp_sparse_JPM2 = sp.SparseMatrix(Sp_JPM2)

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
Sr_sparse_JPM2 = sp.SparseMatrix(Sr_JPM2)

# edge_list_JPM2 = pd.read_csv("../mechanisms/jpm_v2/jpm_v2_EdgeList.csv", comment="!")

# print("creating sparse matrices...")
# Sr_sparse_JPM2, Sp_sparse_JPM2, Svv_sparse_JPM2 = create_sparse(edge_list_JPM2)
# Svv_sparse_JPM2, init_col_del_JPM2 = del_zero_col(Svv_sparse_JPM2)

print("computing dimension of nullspace...")
stoichiometric_invariants_JPM2 = len(Svv_sparse_JPM2.T.nullspace())
Svv_np_JPM2 = np.array(Svv_sparse_JPM2, dtype=float)
rank_Svv_JPM2_np = np.linalg.matrix_rank(Svv_np_JPM2)
dim_leftnull_JPM2 = Svv_sparse_JPM2.shape[0] - rank_Svv_JPM2_np

print("identifying coproduction columns...")
coproduction_cols_JPM2 = create_coproduction(Sr_sparse_JPM2)
print("creating symbolic dictionary...")
symbol_dict_JPM2 = create_symbols(coproduction_cols_JPM2)
print("merging coproduction columns...")
S_merge_JPM2, col_del_JPM2 = merge_coprod(Sr_sparse_JPM2, Sp_sparse_JPM2, symbol_dict_JPM2, coproduction_cols_JPM2, Svv_sparse_JPM2)
S_merge_JPM2[15, 2] = -2*symbol_dict_JPM2[3]
print("performing linear algebra...")
del_r_jpm2 = S_merge_JPM2.shape[1] - Svv_sparse_JPM2.shape[1]

rank_list_jpm2 = linalg_experiment(S_merge_JPM2, num_experiments)
