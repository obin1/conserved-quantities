import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

# no Oxygen tracked
Svv_JPM = sp.Matrix([
               [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [-1, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 1],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, -1, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 1, -1, 2, -1, -1, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
               [-1, 1, 0, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0]])
Svv_sparse_JPM = sp.SparseMatrix(Svv_JPM)

Sp_JPM = sp.Matrix([
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
Sp_sparse_JPM = sp.SparseMatrix(Sp_JPM)

Sr_JPM = sp.Matrix([
               [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 0, -1, 0, -1, -1, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0]])
Sr_sparse_JPM = sp.SparseMatrix(Sr_JPM)

print("computing dimension of nullspace...")
stoichiometric_invariants_JPM = len(Svv_sparse_JPM.T.nullspace())
Svv_np_JPM = np.array(Svv_sparse_JPM, dtype=float)
rank_Svv_JPM_np = np.linalg.matrix_rank(Svv_np_JPM)
dim_leftnull_JPM = Svv_sparse_JPM.shape[0] - rank_Svv_JPM_np

if stoichiometric_invariants_JPM == dim_leftnull_JPM:
    print("identifying coproduction columns...")
    coproduction_cols_JPM = create_coproduction(Sr_sparse_JPM)
    print("creating symbolic dictionary...")
    symbol_dict_JPM = create_symbols(coproduction_cols_JPM)
    print("merging coproduction columns...")
    S_merge_JPM, col_del_JPM = merge_coprod(Sr_sparse_JPM, Sp_sparse_JPM, symbol_dict_JPM, coproduction_cols_JPM, Svv_sparse_JPM)
    S_merge_JPM[15, 2] = -2*symbol_dict_JPM[3]
    print("performing linear algebra...")
    # rank_jpm = S_merge_jpm.rank()
    # dim_null_jpm = S_merge_jpm.shape[0] - rank_jpm
    # del_l_jpm = dim_null_jpm - stoichiometric_invariants_jpm
    del_r_jpm = S_merge_JPM.shape[1] - Svv_sparse_JPM.shape[1]
    # del_c_jpm = -del_r_jpm - del_l_jpm

    rank_list_jpm = linalg_experiment(S_merge_JPM, num_experiments)
    # del_l = S_merge_jpm.shape[0] - rank - stoich_invariants
    # del_c_jpm = -del_r_jpm - del_l_jpm

    # del_l_jpm, del_r_jpm, del_c_jpm = s_linalg(Svv_sparse_jpm, S_merge_jpm, stoichiometric_invariants_jpm)
else:
    print("Error: numpy vs sympy disparity")