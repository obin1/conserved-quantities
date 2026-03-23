import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_gc13_3_0 = pd.read_csv("../mechanisms/geoschem-13.3.0/gckpp_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_gc13_3_0, Sp_sparse_gc13_3_0, Svv_sparse_gc13_3_0 = create_sparse(edge_list_gc13_3_0)
Svv_sparse_gc13_3_0, init_col_del_gc13_3_0 = del_zero_col(Svv_sparse_gc13_3_0)
print("computing dimension of nullspace...")
Svv_gc13_3_0_np = np.array(Svv_sparse_gc13_3_0, dtype=float)
rank_Svv_gc13_3_0_np = np.linalg.matrix_rank(Svv_gc13_3_0_np)
dim_leftnull_gc13_3_0 = Svv_sparse_gc13_3_0.shape[0] - rank_Svv_gc13_3_0_np


print("identifying coproduction columns...")
coproduction_cols_gc13_3_0 = create_coproduction(Sr_sparse_gc13_3_0)
print("creating symbolic dictionary...")
symbol_dict_gc13_3_0 = create_symbols(coproduction_cols_gc13_3_0)
print("merging coproduction columns...")
S_merge_gc13_3_0, col_del_gc13_3_0 = merge_coprod(Sr_sparse_gc13_3_0, Sp_sparse_gc13_3_0, symbol_dict_gc13_3_0, coproduction_cols_gc13_3_0, Svv_sparse_gc13_3_0)
print("performing linear algebra...")

del_r_gc13_3_0 = S_merge_gc13_3_0.shape[1] - Svv_sparse_gc13_3_0.shape[1]
    
rank_list_gc13_3_0 = linalg_experiment(S_merge_gc13_3_0, num_experiments)

A = numeric_sparse_matrix_fast_combined(S_merge_gc13_3_0)
rank_list_gc13_3_0_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])