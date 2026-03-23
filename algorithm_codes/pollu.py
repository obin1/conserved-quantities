import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_pollu = pd.read_csv("../mechanisms/pollu/pollu_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_pollu, Sp_sparse_pollu, Svv_sparse_pollu = create_sparse(edge_list_pollu)
Svv_sparse_pollu, init_col_del_pollu = del_zero_col(Svv_sparse_pollu)
print("computing dimension of nullspace...")
stoichiometric_invariants_pollu = len(Svv_sparse_pollu.T.nullspace())
Svv_pollu_np = np.array(Svv_sparse_pollu, dtype=float)
rank_Svv_pollu_np = np.linalg.matrix_rank(Svv_pollu_np)
dim_leftnull_pollu = Svv_sparse_pollu.shape[0] - rank_Svv_pollu_np

if stoichiometric_invariants_pollu == dim_leftnull_pollu:
    print("identifying coproduction columns...")
    coproduction_cols_pollu = create_coproduction(Sr_sparse_pollu)
    print("creating symbolic dictionary...")
    symbol_dict_pollu = create_symbols(coproduction_cols_pollu)
    print("merging coproduction columns...")
    S_merge_pollu, col_del_pollu = merge_coprod(Sr_sparse_pollu, Sp_sparse_pollu, symbol_dict_pollu, coproduction_cols_pollu, Svv_sparse_pollu)
    print("performing linear algebra...")
    del_r_pollu = S_merge_pollu.shape[1] - Svv_sparse_pollu.shape[1]
    
    rank_list_pollu = linalg_experiment(S_merge_pollu, num_experiments)
    
else:
    print("Error: numpy vs sympy disparity")
