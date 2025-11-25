import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_mcm = pd.read_csv("../mechanisms/mcm_v3.3.1/mcm_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_mcm, Sp_sparse_mcm, Svv_sparse_mcm = create_sparse(edge_list_mcm)
print("computing dimension of nullspace...")
# stoichiometric_invariants_gc14_6_3 = len(Svv_sparse_gc14_6_3.T.nullspace())
Svv_mcm_np = np.array(Svv_sparse_mcm, dtype=float)
rank_Svv_mcm_np = np.linalg.matrix_rank(Svv_mcm_np)
dim_leftnull_mcm = Svv_sparse_mcm.shape[0] - rank_Svv_mcm_np


print("identifying coproduction columns...")
coproduction_cols_mcm = create_coproduction(Sr_sparse_mcm)
print("creating symbolic dictionary...")
symbol_dict_mcm = create_symbols(coproduction_cols_mcm)
print("merging coproduction columns...")
S_merge_mcm, col_del_mcm = merge_coprod(Sr_sparse_mcm, Sp_sparse_mcm, symbol_dict_mcm, coproduction_cols_mcm, Svv_sparse_mcm)
print("performing linear algebra...")

del_r_mcm = S_merge_mcm.shape[1] - Svv_sparse_mcm.shape[1]
    
rank_list_mcm = linalg_experiment(S_merge_mcm, num_experiments)

