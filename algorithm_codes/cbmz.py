import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cbmz = pd.read_csv("../mechanisms/cbmz/CBMZ_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cbmz, Sp_sparse_cbmz, Svv_sparse_cbmz = create_sparse(edge_list_cbmz)
Svv_sparse_cbmz, init_col_del_cbmz = del_zero_col(Svv_sparse_cbmz)
print("computing dimension of nullspace...")
stoichiometric_invariants_cbmz = len(Svv_sparse_cbmz.T.nullspace())
Svv_cbmz_np = np.array(Svv_sparse_cbmz, dtype=float)
rank_Svv_cbmz_np = np.linalg.matrix_rank(Svv_cbmz_np)
dim_leftnull_cbmz = Svv_sparse_cbmz.shape[0] - rank_Svv_cbmz_np

print("identifying coproduction columns...")
coproduction_cols_cbmz = create_coproduction(Sr_sparse_cbmz)
print("creating symbolic dictionary...")
symbol_dict_cbmz = create_symbols(coproduction_cols_cbmz)
print("merging coproduction columns...")
S_merge_cbmz, col_del_cbmz = merge_coprod(Sr_sparse_cbmz, Sp_sparse_cbmz, symbol_dict_cbmz, coproduction_cols_cbmz, Svv_sparse_cbmz)
print("performing linear algebra...")

del_r_cbmz = S_merge_cbmz.shape[1] - Svv_sparse_cbmz.shape[1]

rank_list_cbmz = linalg_experiment(S_merge_cbmz, num_experiments)