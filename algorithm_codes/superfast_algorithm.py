import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

edge_list_superfast = pd.read_csv("../mechanisms/superfast/superfast_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_superfast, Sp_sparse_superfast, Svv_sparse_superfast = create_sparse(edge_list_superfast)
print("computing dimension of nullspace...")
stoichiometric_invariants_superfast = len(Svv_sparse_superfast.T.nullspace())
Svv_superfast_np = np.array(Svv_sparse_superfast, dtype=float)
rank_Svv_superfast_np = np.linalg.matrix_rank(Svv_superfast_np)
dim_leftnull_superfast = Svv_sparse_superfast.shape[0] - rank_Svv_superfast_np

if stoichiometric_invariants_superfast == dim_leftnull_superfast:
    print("identifying coproduction columns...")
    coproduction_cols_superfast = create_coproduction(Sr_sparse_superfast)
    print("creating symbolic dictionary...")
    symbol_dict_superfast = create_symbols(coproduction_cols_superfast)
    print("merging coproduction columns...")
    S_merge_superfast = merge_coprod(Sr_sparse_superfast, Sp_sparse_superfast, symbol_dict_superfast, coproduction_cols_superfast, Svv_sparse_superfast)
    print("performing linear algebra...")
    del_l_superfast, del_r_superfast, del_c_superfast = s_linalg(Svv_sparse_superfast, S_merge_superfast, stoichiometric_invariants_superfast)
else:
    print("Error: numpy vs sympy disparity")
