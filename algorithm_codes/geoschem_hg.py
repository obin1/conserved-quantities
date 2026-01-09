import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_gchg= pd.read_csv("../mechanisms/geoschem-hg/Hg_EdgeList.csv", comment="!")
print("creating sparse matrices...")
Sr_sparse_gchg, Sp_sparse_gchg, Svv_sparse_gchg= create_sparse(edge_list_gchg)
Svv_sparse_gchg, init_col_del_gchg = del_zero_col(Svv_sparse_gchg)
print("computing dimension of nullspace...")
# stoichiometric_invariants_gchg= len(Svv_sparse_gchg.T.nullspace())
Svv_gchg_np = np.array(Svv_sparse_gchg, dtype=float)
rank_Svv_gchg_np = np.linalg.matrix_rank(Svv_gchg_np)
dim_leftnull_gchg= Svv_sparse_gchg.shape[0] - rank_Svv_gchg_np


print("identifying coproduction columns...")
coproduction_cols_gchg= create_coproduction(Sr_sparse_gchg)
print("creating symbolic dictionary...")
symbol_dict_gchg= create_symbols(coproduction_cols_gchg)
print("merging coproduction columns...")
S_merge_gchg, col_del_gchg= merge_coprod(Sr_sparse_gchg, Sp_sparse_gchg, symbol_dict_gchg, coproduction_cols_gchg, Svv_sparse_gchg)
print("performing linear algebra...")

del_r_gchg= S_merge_gchg.shape[1] - Svv_sparse_gchg.shape[1]
    
rank_list_gchg= linalg_experiment(S_merge_gchg, num_experiments)

