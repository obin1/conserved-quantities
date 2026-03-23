import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment
from scipy.linalg import null_space

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_gc = pd.read_csv("../mechanisms/geos-chem-v14/gckpp_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_gc, Sp_sparse_gc, Svv_sparse_gc = create_sparse(edge_list_gc)
Svv_sparse_gc, init_col_del_gc = del_zero_col(Svv_sparse_gc)
print("computing dimension of nullspace...")
gc_null = Svv_sparse_gc.T.nullspace()
stoichiometric_invariants_gc = len(gc_null)

Svv_gc_np = np.array(Svv_sparse_gc, dtype=float)
rank_Svv_gc_np = np.linalg.matrix_rank(Svv_gc_np)
dim_leftnull_gc = Svv_gc_np.shape[0] - rank_Svv_gc_np
gc_null_np = null_space(Svv_gc_np.T)

print("identifying coproduction columns...")
coproduction_cols_gc = create_coproduction(Sr_sparse_gc)
print("creating symbolic dictionary...")
symbol_dict_gc = create_symbols(coproduction_cols_gc)
print("merging coproduction columns...")
S_merge_gc, col_del_gc = merge_coprod(Sr_sparse_gc, Sp_sparse_gc, symbol_dict_gc, coproduction_cols_gc, Svv_sparse_gc)
print("performing linear algebra...")

del_r_gc = S_merge_gc.shape[1] - Svv_sparse_gc.shape[1]

rank_list_gc = linalg_experiment(S_merge_gc, num_experiments)
