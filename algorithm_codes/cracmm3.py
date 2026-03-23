# WITH DELTA_C, DELTA_N, DELTA_SI
import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cracmm3 = pd.read_csv("../mechanisms/cracmm/cracmm3_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cracmm3, Sp_sparse_cracmm3, Svv_sparse_cracmm3 = create_sparse(edge_list_cracmm3)
Svv_sparse_cracmm3, init_col_del_cracmm3 = del_zero_col(Svv_sparse_cracmm3)
print("computing dimension of nullspace...")
stoichiometric_invariants_cracmm3 = len(Svv_sparse_cracmm3.T.nullspace())
Svv_cracmm3_np = np.array(Svv_sparse_cracmm3, dtype=float)
rank_Svv_cracmm3_np = np.linalg.matrix_rank(Svv_cracmm3_np)
dim_leftnull_cracmm3 = Svv_sparse_cracmm3.shape[0] - rank_Svv_cracmm3_np

print("identifying coproduction columns...")
coproduction_cols_cracmm3 = create_coproduction(Sr_sparse_cracmm3)
print("creating symbolic dictionary...")
symbol_dict_cracmm3 = create_symbols(coproduction_cols_cracmm3)
print("merging coproduction columns...")
S_merge_cracmm3, col_del_cracmm3 = merge_coprod(Sr_sparse_cracmm3, Sp_sparse_cracmm3, symbol_dict_cracmm3, coproduction_cols_cracmm3, Svv_sparse_cracmm3)
print("performing linear algebra...")
del_r_cracmm3 = S_merge_cracmm3.shape[1] - Svv_sparse_cracmm3.shape[1]
    
rank_list_cracmm3 = linalg_experiment(S_merge_cracmm3, num_experiments)


#%%
# WITHOUT DELTA_C, DELTA_N, DELTA_SI

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cracmm3e = pd.read_csv("../mechanisms/cracmm/cracmm3_EdgeList_eliminated.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cracmm3e, Sp_sparse_cracmm3e, Svv_sparse_cracmm3e = create_sparse(edge_list_cracmm3e)
Svv_sparse_cracmm3e, init_col_del_cracmm3e = del_zero_col(Svv_sparse_cracmm3e)
print("computing dimension of nullspace...")
stoichiometric_invariants_cracmm3e = len(Svv_sparse_cracmm3e.T.nullspace())
Svv_cracmm3e_np = np.array(Svv_sparse_cracmm3e, dtype=float)
rank_Svv_cracmm3e_np = np.linalg.matrix_rank(Svv_cracmm3e_np)
dim_leftnull_cracmm3e = Svv_sparse_cracmm3e.shape[0] - rank_Svv_cracmm3e_np

print("identifying coproduction columns...")
coproduction_cols_cracmm3e = create_coproduction(Sr_sparse_cracmm3e)
print("creating symbolic dictionary...")
symbol_dict_cracmm3e = create_symbols(coproduction_cols_cracmm3e)
print("merging coproduction columns...")
S_merge_cracmm3e, col_del_cracmm3e = merge_coprod(Sr_sparse_cracmm3e, Sp_sparse_cracmm3e, symbol_dict_cracmm3e, coproduction_cols_cracmm3e, Svv_sparse_cracmm3e)
print("performing linear algebra...")
del_r_cracmm3e = S_merge_cracmm3e.shape[1] - Svv_sparse_cracmm3e.shape[1]
    
rank_list_cracmm3e = linalg_experiment(S_merge_cracmm3e, num_experiments)
