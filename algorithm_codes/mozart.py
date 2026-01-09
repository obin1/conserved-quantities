import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_mozart = pd.read_csv("../mechanisms/mozart/MOZART_4_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_mozart, Sp_sparse_mozart, Svv_sparse_mozart = create_sparse(edge_list_mozart)
Svv_sparse_mozart, init_col_del_mozart = del_zero_col(Svv_sparse_mozart)
print("computing dimension of nullspace...")
stoichiometric_invariants_mozart = len(Svv_sparse_mozart.T.nullspace())
Svv_mozart_np = np.array(Svv_sparse_mozart, dtype=float)
rank_Svv_mozart_np = np.linalg.matrix_rank(Svv_mozart_np)
dim_leftnull_mozart = Svv_sparse_mozart.shape[0] - rank_Svv_mozart_np

if stoichiometric_invariants_mozart == dim_leftnull_mozart:
    print("identifying coproduction columns...")
    coproduction_cols_mozart = create_coproduction(Sr_sparse_mozart)
    print("creating symbolic dictionary...")
    symbol_dict_mozart = create_symbols(coproduction_cols_mozart)
    print("merging coproduction columns...")
    S_merge_mozart, col_del_mozart = merge_coprod(
        Sr_sparse_mozart,
        Sp_sparse_mozart,
        symbol_dict_mozart,
        coproduction_cols_mozart,
        Svv_sparse_mozart
    )
    print("performing linear algebra...")
    # rank_mozart = S_merge_mozart.rank()
    # dim_null_mozart = S_merge_mozart.shape[0] - rank_mozart
    # del_l_mozart = dim_null_mozart - stoichiometric_invariants_mozart
    del_r_mozart = S_merge_mozart.shape[1] - Svv_sparse_mozart.shape[1]
    # del_c_mozart = -del_r_mozart - del_l_mozart
    
    rank_list_mozart = linalg_experiment(S_merge_mozart, num_experiments)
    # del_l_mozart = S_merge_mozart.shape[0] - rank - stoich_invariants
    # del_c_mozart = -del_r_mozart - del_l_mozart

    # del_l_mozart, del_r_mozart, del_c_mozart = s_linalg(
    #     Svv_sparse_mozart,
    #     S_merge_mozart,
    #     stoichiometric_invariants_mozart
    # )
else:
    print("Error: numpy vs sympy disparity")
