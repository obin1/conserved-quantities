import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

pd.set_option('display.max_columns', None)
edge_list_logan = pd.read_csv("../mechanisms/logan/log81_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_logan, Sp_sparse_logan, Svv_sparse_logan = create_sparse(edge_list_logan)
Svv_sparse_logan, init_col_del_logan = del_zero_col(Svv_sparse_logan)
print("computing dimension of nullspace...")
# stoichiometric_invariants_logan = len(Svv_sparse_logan.T.nullspace())
Svv_logan_np = np.array(Svv_sparse_logan, dtype=float)
rank_Svv_logan_np = np.linalg.matrix_rank(Svv_logan_np)
dim_leftnull_logan = Svv_sparse_logan.shape[0] - rank_Svv_logan_np


print("identifying coproduction columns...")
coproduction_cols_logan = create_coproduction(Sr_sparse_logan)
print("creating symbolic dictionary...")
symbol_dict_logan = create_symbols(coproduction_cols_logan)
print("merging coproduction columns...")
S_merge_logan, col_del_logan = merge_coprod(Sr_sparse_logan, Sp_sparse_logan, symbol_dict_logan, coproduction_cols_logan, Svv_sparse_logan)
print("performing linear algebra...")

del_r_logan = S_merge_logan.shape[1] - Svv_sparse_logan.shape[1]
    
rank_list_logan = linalg_experiment(S_merge_logan, num_experiments)

