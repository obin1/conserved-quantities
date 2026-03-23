import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_mecca = pd.read_csv("../mechanisms/mecca_v4.6.0/gas_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_mecca, Sp_sparse_mecca, Svv_sparse_mecca = create_sparse(edge_list_mecca)
Svv_sparse_mecca, init_col_del_mecca = del_zero_col(Svv_sparse_mecca)
print("computing dimension of nullspace...")
Svv_mecca_np = np.array(Svv_sparse_mecca, dtype=float)
rank_Svv_mecca_np = np.linalg.matrix_rank(Svv_mecca_np)
dim_leftnull_mecca = Svv_sparse_mecca.shape[0] - rank_Svv_mecca_np


print("identifying coproduction columns...")
coproduction_cols_mecca = create_coproduction(Sr_sparse_mecca)
print("creating symbolic dictionary...")
symbol_dict_mecca = create_symbols(coproduction_cols_mecca)
print("merging coproduction columns...")
S_merge_mecca, col_del_mecca = merge_coprod(Sr_sparse_mecca, Sp_sparse_mecca, symbol_dict_mecca, coproduction_cols_mecca, Svv_sparse_mecca)
print("performing linear algebra...")
del_r_mecca = S_merge_mecca.shape[1] - Svv_sparse_mecca.shape[1]
    
rank_list_mecca = linalg_experiment(S_merge_mecca, num_experiments)

