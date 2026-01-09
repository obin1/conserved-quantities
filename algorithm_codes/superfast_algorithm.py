import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_superfast = pd.read_csv("../mechanisms/superfast/superfast_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_superfast, Sp_sparse_superfast, Svv_sparse_superfast = create_sparse(edge_list_superfast)
Svv_sparse_superfast, init_col_del_superfast = del_zero_col(Svv_sparse_superfast)
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
    S_merge_superfast, col_del_supefast = merge_coprod(Sr_sparse_superfast, Sp_sparse_superfast, symbol_dict_superfast, coproduction_cols_superfast, Svv_sparse_superfast)
    print("performing linear algebra...")
    # rank_superfast = S_merge_superfast.rank()
    # dim_null_superfast = S_merge_superfast.shape[0] - rank_superfast
    # del_l_superfast = dim_null_superfast - stoichiometric_invariants_superfast
    del_r_superfast = S_merge_superfast.shape[1] - Svv_sparse_superfast.shape[1]
    # del_c_superfast = -del_r_superfast - del_l_superfast

    rank_list_superfast = linalg_experiment(S_merge_superfast, num_experiments)
    # del_l = S_merge_superfast.shape[0] - rank - stoich_invariants
    # del_c_superfast = -del_r_superfast - del_l_superfast

    # del_l_superfast, del_r_superfast, del_c_superfast = s_linalg(Svv_sparse_superfast, S_merge_superfast, stoichiometric_invariants_superfast)
else:
    print("Error: numpy vs sympy disparity")
