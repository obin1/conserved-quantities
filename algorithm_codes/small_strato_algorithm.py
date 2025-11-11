import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

edge_list_small_strato = pd.read_csv("../mechanisms/small_strato/small_strato_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_small_strato, Sp_sparse_small_strato, Svv_sparse_small_strato = create_sparse(edge_list_small_strato)
print("computing dimension of nullspace...")
stoichiometric_invariants_small_strato = len(Svv_sparse_small_strato.T.nullspace())
Svv_small_strato_np = np.array(Svv_sparse_small_strato, dtype=float)
rank_Svv_small_strato_np = np.linalg.matrix_rank(Svv_small_strato_np)
dim_leftnull_small_strato = Svv_sparse_small_strato.shape[0] - rank_Svv_small_strato_np

if stoichiometric_invariants_small_strato == dim_leftnull_small_strato:
    print("identifying coproduction columns...")
    coproduction_cols_small_strato = create_coproduction(Sr_sparse_small_strato)
    print("creating symbolic dictionary...")
    symbol_dict_small_strato = create_symbols(coproduction_cols_small_strato)
    print("merging coproduction columns...")
    S_merge_small_strato = merge_coprod(Sr_sparse_small_strato, Sp_sparse_small_strato, symbol_dict_small_strato, coproduction_cols_small_strato, Svv_sparse_small_strato)
    print("performing linear algebra...")
    del_l_small_strato, del_r_small_strato, del_c_small_strato = s_linalg(Svv_sparse_small_strato, S_merge_small_strato, stoichiometric_invariants_small_strato)
else:
    print("Error: numpy vs sympy disparity")
