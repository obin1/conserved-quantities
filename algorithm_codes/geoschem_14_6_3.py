import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_gc14_6_3 = pd.read_csv("../mechanisms/geoschem-14.6.3/gckpp_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_gc14_6_3, Sp_sparse_gc14_6_3, Svv_sparse_gc14_6_3 = create_sparse(edge_list_gc14_6_3)
Svv_sparse_gc14_6_3, init_col_del_gc14_6_3 = del_zero_col(Svv_sparse_gc14_6_3)
print("computing dimension of nullspace...")
# stoichiometric_invariants_gc14_6_3 = len(Svv_sparse_gc14_6_3.T.nullspace())
Svv_gc14_6_3_np = np.array(Svv_sparse_gc14_6_3, dtype=float)
rank_Svv_gc14_6_3_np = np.linalg.matrix_rank(Svv_gc14_6_3_np)
dim_leftnull_gc14_6_3 = Svv_sparse_gc14_6_3.shape[0] - rank_Svv_gc14_6_3_np


print("identifying coproduction columns...")
coproduction_cols_gc14_6_3 = create_coproduction(Sr_sparse_gc14_6_3)
print("creating symbolic dictionary...")
symbol_dict_gc14_6_3 = create_symbols(coproduction_cols_gc14_6_3)
print("merging coproduction columns...")
S_merge_gc14_6_3, col_del_gc14_6_3 = merge_coprod(Sr_sparse_gc14_6_3, Sp_sparse_gc14_6_3, symbol_dict_gc14_6_3, coproduction_cols_gc14_6_3, Svv_sparse_gc14_6_3)
print("performing linear algebra...")

del_r_gc14_6_3 = S_merge_gc14_6_3.shape[1] - Svv_sparse_gc14_6_3.shape[1]
    
rank_list_gc14_6_3 = linalg_experiment(S_merge_gc14_6_3, num_experiments)

