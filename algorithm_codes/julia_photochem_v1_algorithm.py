import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_JPM = pd.read_csv("../mechanisms/jpm_v1/jpm_v1_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_JPM, Sp_sparse_JPM, Svv_sparse_JPM = create_sparse(edge_list_JPM)
Svv_sparse_JPM, init_col_del_JPM = del_zero_col(Svv_sparse_JPM)

print("computing dimension of nullspace...")
stoichiometric_invariants_JPM = len(Svv_sparse_JPM.T.nullspace())
Svv_np_JPM = np.array(Svv_sparse_JPM, dtype=float)
rank_Svv_JPM_np = np.linalg.matrix_rank(Svv_np_JPM)
dim_leftnull_JPM = Svv_sparse_JPM.shape[0] - rank_Svv_JPM_np

print("identifying coproduction columns...")
coproduction_cols_JPM = create_coproduction(Sr_sparse_JPM)
print("creating symbolic dictionary...")
symbol_dict_JPM = create_symbols(coproduction_cols_JPM)
print("merging coproduction columns...")
S_merge_JPM, col_del_JPM = merge_coprod(Sr_sparse_JPM, Sp_sparse_JPM, symbol_dict_JPM, coproduction_cols_JPM, Svv_sparse_JPM)
# S_merge_JPM[15, 2] = -2*symbol_dict_JPM[3]
print("performing linear algebra...")
del_r_jpm = S_merge_JPM.shape[1] - Svv_sparse_JPM.shape[1]

rank_list_jpm = linalg_experiment(S_merge_JPM, num_experiments)
