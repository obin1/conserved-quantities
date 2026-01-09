import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_mozart_t1 = pd.read_csv("../mechanisms/mozart-t1/MOZART_T1_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_mozart_t1, Sp_sparse_mozart_t1, Svv_sparse_mozart_t1 = create_sparse(edge_list_mozart_t1)
Svv_sparse_mozart_t1, init_col_del_mozart_t1 = del_zero_col(Svv_sparse_mozart_t1)
print("computing dimension of nullspace...")
stoichiometric_invariants_mozart_t1 = len(Svv_sparse_mozart_t1.T.nullspace())
Svv_mozart_t1_np = np.array(Svv_sparse_mozart_t1, dtype=float)
rank_Svv_mozart_t1_np = np.linalg.matrix_rank(Svv_mozart_t1_np)
dim_leftnull_mozart_t1 = Svv_sparse_mozart_t1.shape[0] - rank_Svv_mozart_t1_np

if stoichiometric_invariants_mozart_t1 == dim_leftnull_mozart_t1:
    print("identifying coproduction columns...")
    coproduction_cols_mozart_t1 = create_coproduction(Sr_sparse_mozart_t1)
    print("creating symbolic dictionary...")
    symbol_dict_mozart_t1 = create_symbols(coproduction_cols_mozart_t1)
    print("merging coproduction columns...")
    S_merge_mozart_t1, col_del_mozart_t1 = merge_coprod(
        Sr_sparse_mozart_t1,
        Sp_sparse_mozart_t1,
        symbol_dict_mozart_t1,
        coproduction_cols_mozart_t1,
        Svv_sparse_mozart_t1
    )
    print("performing linear algebra...")
    # rank_mozart_t1 = S_merge_mozart_t1.rank()
    # dim_null_mozart_t1 = S_merge_mozart_t1.shape[0] - rank_mozart_t1
    # del_l_mozart_t1 = dim_null_mozart_t1 - stoichiometric_invariants_mozart_t1
    del_r_mozart_t1 = S_merge_mozart_t1.shape[1] - Svv_sparse_mozart_t1.shape[1]
    # del_c_mozart_t1 = -del_r_mozart_t1 - del_l_mozart_t1
    
    rank_list_mozart_t1 = linalg_experiment(S_merge_mozart_t1, num_experiments)
    # del_l_mozart_t1 = S_merge_mozart_t1.shape[0] - rank - stoich_invariants
    # del_c_mozart_t1 = -del_r_mozart_t1 - del_l_mozart_t1

    # del_l_mozart_t1, del_r_mozart_t1, del_c_mozart_t1 = s_linalg(
    #     Svv_sparse_mozart_t1,
    #     S_merge_mozart_t1,
    #     stoichiometric_invariants_mozart_t1
    # )
else:
    print("Error: numpy vs sympy disparity")
