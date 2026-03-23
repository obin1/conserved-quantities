import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_jam = pd.read_csv("../mechanisms/jam/jam_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_jam, Sp_sparse_jam, Svv_sparse_jam = create_sparse(edge_list_jam)
Svv_sparse_jam, init_col_del_jam = del_zero_col(Svv_sparse_jam)
print("computing dimension of nullspace...")
Svv_jam_np = np.array(Svv_sparse_jam, dtype=float)
rank_Svv_jam_np = np.linalg.matrix_rank(Svv_jam_np)
dim_leftnull_jam = Svv_sparse_jam.shape[0] - rank_Svv_jam_np


print("identifying coproduction columns...")
coproduction_cols_jam = create_coproduction(Sr_sparse_jam)
print("creating symbolic dictionary...")
symbol_dict_jam = create_symbols(coproduction_cols_jam)
print("merging coproduction columns...")
S_merge_jam, col_del_jam = merge_coprod(Sr_sparse_jam, Sp_sparse_jam, symbol_dict_jam, coproduction_cols_jam, Svv_sparse_jam)
print("performing linear algebra...")

del_r_jam = S_merge_jam.shape[1] - Svv_sparse_jam.shape[1]
    
rank_list_jam = linalg_experiment(S_merge_jam, num_experiments)

