import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 1

edge_list_gc12_8_0 = pd.read_csv("../mechanisms/geoschem-12.8.0/gckpp_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_gc12_8_0, Sp_sparse_gc12_8_0, Svv_sparse_gc12_8_0 = create_sparse(edge_list_gc12_8_0)
print("computing dimension of nullspace...")
# stoichiometric_invariants_gc12_8_0 = len(Svv_sparse_gc12_8_0.T.nullspace())
Svv_gc12_8_0_np = np.array(Svv_sparse_gc12_8_0, dtype=float)
rank_Svv_gc12_8_0_np = np.linalg.matrix_rank(Svv_gc12_8_0_np)
dim_leftnull_gc12_8_0 = Svv_sparse_gc12_8_0.shape[0] - rank_Svv_gc12_8_0_np


print("identifying coproduction columns...")
coproduction_cols_gc12_8_0 = create_coproduction(Sr_sparse_gc12_8_0)
print("creating symbolic dictionary...")
symbol_dict_gc12_8_0 = create_symbols(coproduction_cols_gc12_8_0)
print("merging coproduction columns...")
S_merge_gc12_8_0, col_del_gc12_8_0 = merge_coprod(Sr_sparse_gc12_8_0, Sp_sparse_gc12_8_0, symbol_dict_gc12_8_0, coproduction_cols_gc12_8_0, Svv_sparse_gc12_8_0)
print("performing linear algebra...")

del_r_gc12_8_0 = S_merge_gc12_8_0.shape[1] - Svv_sparse_gc12_8_0.shape[1]
    
rank_list_gc12_8_0 = linalg_experiment(S_merge_gc12_8_0, num_experiments)

