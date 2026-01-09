import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, create_coproduction_2, linalg_experiment_fast, numeric_sparse_matrix_fast_combined

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_mcm = pd.read_csv("../mechanisms/mcm_v3.3.1/mcm_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_mcm, Sp_sparse_mcm, Svv_sparse_mcm = create_sparse(edge_list_mcm)
Svv_sparse_mcm, init_col_del_mcm = del_zero_col(Svv_sparse_mcm)
print("computing dimension of nullspace...")
Svv_mcm_np = np.array(Svv_sparse_mcm, dtype=float)
rank_Svv_mcm_np = np.linalg.matrix_rank(Svv_mcm_np)
dim_leftnull_mcm = Svv_sparse_mcm.shape[0] - rank_Svv_mcm_np


print("identifying coproduction columns...")
# coproduction_cols_mcm = create_coproduction(Sr_sparse_mcm)
coproduction_cols_mcm = create_coproduction_2(Sr_sparse_mcm)
print("creating symbolic dictionary...")
symbol_dict_mcm = create_symbols(coproduction_cols_mcm)
print("merging coproduction columns...")
S_merge_mcm, col_del_mcm = merge_coprod(Sr_sparse_mcm, Sp_sparse_mcm, symbol_dict_mcm, coproduction_cols_mcm, Svv_sparse_mcm)
print("performing linear algebra...")

del_r_mcm = S_merge_mcm.shape[1] - Svv_sparse_mcm.shape[1]
    
# rank_list_mcm = linalg_experiment(S_merge_mcm, num_experiments)
# rank_list_mcm = linalg_experiment_fast(S_merge_mcm, num_experiments, proj_dim=300)

A = numeric_sparse_matrix_fast_combined(S_merge_mcm)
rank_list_mcm = linalg_experiment_fast(A, num_experiments, A.shape[0])
# del_l = S_merge_mcm.shape[0] - rank - stoich_invariants
# del_c_mcm = -del_r_mcm - del_l_mcm

