import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_racm = pd.read_csv("../mechanisms/racm/RACM_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_racm, Sp_sparse_racm, Svv_sparse_racm = create_sparse(edge_list_racm)
Svv_sparse_racm, init_col_del_racm = del_zero_col(Svv_sparse_racm)
print("computing dimension of nullspace...")
stoichiometric_invariants_racm = len(Svv_sparse_racm.T.nullspace())
Svv_racm_np = np.array(Svv_sparse_racm, dtype=float)
rank_Svv_racm_np = np.linalg.matrix_rank(Svv_racm_np)
dim_leftnull_racm = Svv_sparse_racm.shape[0] - rank_Svv_racm_np

if stoichiometric_invariants_racm == dim_leftnull_racm:
    print("identifying coproduction columns...")
    coproduction_cols_racm = create_coproduction(Sr_sparse_racm)
    print("creating symbolic dictionary...")
    symbol_dict_racm = create_symbols(coproduction_cols_racm)
    print("merging coproduction columns...")
    S_merge_racm, col_del_racm = merge_coprod(Sr_sparse_racm, Sp_sparse_racm, symbol_dict_racm, coproduction_cols_racm, Svv_sparse_racm)
    print("performing linear algebra...")
    del_r_racm = S_merge_racm.shape[1] - Svv_sparse_racm.shape[1]
    
    rank_list_racm = linalg_experiment(S_merge_racm, num_experiments)

else:
    print("Error: numpy vs sympy disparity")
