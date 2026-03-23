import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_saprc99 = pd.read_csv("../mechanisms/saprc99/saprc99_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_saprc99, Sp_sparse_saprc99, Svv_sparse_saprc99 = create_sparse(edge_list_saprc99)
Svv_sparse_saprc99, init_col_del_saprc99 = del_zero_col(Svv_sparse_saprc99)
print("computing dimension of nullspace...")
stoichiometric_invariants_saprc99 = len(Svv_sparse_saprc99.T.nullspace())
Svv_saprc99_np = np.array(Svv_sparse_saprc99, dtype=float)
rank_Svv_saprc99_np = np.linalg.matrix_rank(Svv_saprc99_np)
dim_leftnull_saprc99 = Svv_sparse_saprc99.shape[0] - rank_Svv_saprc99_np

if stoichiometric_invariants_saprc99 == dim_leftnull_saprc99:
    print("identifying coproduction columns...")
    coproduction_cols_saprc99 = create_coproduction(Sr_sparse_saprc99)
    print("creating symbolic dictionary...")
    symbol_dict_saprc99 = create_symbols(coproduction_cols_saprc99)
    print("merging coproduction columns...")
    S_merge_saprc99, col_del_saprc99 = merge_coprod(Sr_sparse_saprc99, Sp_sparse_saprc99, symbol_dict_saprc99, coproduction_cols_saprc99, Svv_sparse_saprc99)
    print("performing linear algebra...")
    del_r_saprc99 = S_merge_saprc99.shape[1] - Svv_sparse_saprc99.shape[1]

    rank_list_saprc99 = linalg_experiment(S_merge_saprc99, num_experiments)
    
else:
    print("Error: numpy vs sympy disparity")
