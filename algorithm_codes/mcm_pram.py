import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, create_coproduction_2, linalg_experiment_fast, numeric_sparse_matrix_fast_combined

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_pram = pd.read_csv("../mechanisms/mcm-pram/mcm_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_pram, Sp_sparse_pram, Svv_sparse_pram = create_sparse(edge_list_pram)
Svv_sparse_pram, init_col_del_pram = del_zero_col(Svv_sparse_pram)
print("computing dimension of nullspace...")
Svv_pram_np = np.array(Svv_sparse_pram, dtype=float)
rank_Svv_pram_np = np.linalg.matrix_rank(Svv_pram_np)
dim_leftnull_pram = Svv_sparse_pram.shape[0] - rank_Svv_pram_np
print("rank_Svv_pram_np:", rank_Svv_pram_np)

print("identifying coproduction columns...")
coproduction_cols_pram = create_coproduction_2(Sr_sparse_pram)
print("creating symbolic dictionary...")
symbol_dict_pram = create_symbols(coproduction_cols_pram)
print("merging coproduction columns...")
S_merge_pram, col_del_pram = merge_coprod(Sr_sparse_pram, Sp_sparse_pram, symbol_dict_pram, coproduction_cols_pram, Svv_sparse_pram)
print("performing linear algebra...")
del_r_pram = S_merge_pram.shape[1] - Svv_sparse_pram.shape[1]
    
A = numeric_sparse_matrix_fast_combined(S_merge_pram)
rank_list_pram = linalg_experiment_fast(A, num_experiments, A.shape[0])
print("rank_list_pram:", rank_list_pram)

