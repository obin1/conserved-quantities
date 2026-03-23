import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment


SEED = 42
np.random.seed(SEED)

num_experiments = 10

Svv_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [-2, -1,  0,  0, -1,  1, -1,  1,  1,  0], # O2
    [ 1, -2, -2,  0,  0, -1,  0, -2,  0,  0], # PhCH2O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 1,  2,  0,  0,  1, -1,  1,  0, -2,  0], # HO2
    [ 0,  2,  1, -1,  1,  0,  1,  0,  0,  0], # PhCHO
    [ 0,  0,  1,  0, -1,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0, -1,  0,  1,  0,  0,  2], # OH
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  1, -1,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1], # H2O2
])
Svv_sparse_d = sp.SparseMatrix(Svv_d)

Sr_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [-2, -1,  0,  0, -1,  0, -1,  0,  0,  0], # O2
    [ 0, -2, -2,  0,  0, -1,  0, -2,  0,  0], # PhCH2O2
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 0,  0,  0,  0,  0, -1,  0,  0, -2,  0], # HO2
    [ 0,  0,  0, -1,  0,  0,  0,  0,  0,  0], # PhCHO
    [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # OH
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # H2O2
])
Sr_sparse_d = sp.SparseMatrix(Sr_d)

Sp_d = sp.Matrix([
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [ 0,  0,  0,  0,  0,  1,  0,  1,  1,  0], # O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 1,  2,  0,  0,  1,  0,  1,  0,  0,  0], # HO2
    [ 0,  2,  1,  0,  1,  0,  1,  0,  0,  0], # PhCHO
    [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  2], # OH
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # H2O2
])
Sp_sparse_d = sp.SparseMatrix(Sp_d)

print("computing dimension of nullspace...")
stoichiometric_invariants_d = len(Svv_sparse_d.T.nullspace())
Svv_d_np = np.array(Svv_sparse_d, dtype=float)
rank_Svv_d_np = np.linalg.matrix_rank(Svv_d_np)
dim_leftnull_d = Svv_sparse_d.shape[0] - rank_Svv_d_np

print("identifying coproduction columns...")
coproduction_cols_d = create_coproduction(Sr_sparse_d)
print("creating symbolic dictionary...")
symbol_dict_d = create_symbols(coproduction_cols_d)
print("merging coproduction columns...")
S_merge_d, col_del_d = merge_coprod(Sr_sparse_d, Sp_sparse_d, symbol_dict_d, coproduction_cols_d, Svv_sparse_d)
print("performing linear algebra...")
del_r_d = S_merge_d.shape[1] - Svv_sparse_d.shape[1]

rank_list_d = linalg_experiment(S_merge_d, num_experiments)