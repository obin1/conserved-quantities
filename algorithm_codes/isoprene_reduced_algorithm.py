import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_iso_red = pd.read_csv("../mechanisms/isoprene_reduced_plus_v5/isoprene_reduced_plus_v5_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_iso_red, Sp_sparse_iso_red, Svv_sparse_iso_red = create_sparse(edge_list_iso_red)
Svv_sparse_iso_red, init_col_del_iso_red = del_zero_col(Svv_sparse_iso_red)
print("computing dimension of nullspace...")
stoichiometric_invariants_iso_red = len(Svv_sparse_iso_red.T.nullspace())
Svv_iso_red_np = np.array(Svv_sparse_iso_red, dtype=float)
rank_Svv_iso_red_np = np.linalg.matrix_rank(Svv_iso_red_np)
dim_leftnull_iso_red = Svv_sparse_iso_red.shape[0] - rank_Svv_iso_red_np

if stoichiometric_invariants_iso_red == dim_leftnull_iso_red:
    print("identifying coproduction columns...")
    coproduction_cols_iso_red = create_coproduction(Sr_sparse_iso_red)
    print("creating symbolic dictionary...")
    symbol_dict_iso_red = create_symbols(coproduction_cols_iso_red)
    print("merging coproduction columns...")
    S_merge_iso_red, col_del_iso_red = merge_coprod(Sr_sparse_iso_red, Sp_sparse_iso_red, symbol_dict_iso_red, coproduction_cols_iso_red, Svv_sparse_iso_red)
    print("performing linear algebra...")
    # rank_iso_red = S_merge_iso_red.rank()
    # dim_null_iso_red = S_merge_iso_red.shape[0] - rank_iso_red
    # del_l_iso_red = dim_null_iso_red - stoichiometric_invariants_iso_red
    del_r_iso_red = S_merge_iso_red.shape[1] - Svv_sparse_iso_red.shape[1]
    # del_c_iso_red = -del_r_iso_red - del_l_iso_red

    rank_list_iso_red = linalg_experiment(S_merge_iso_red, num_experiments)
    # del_l = S_merge_iso_red.shape[0] - rank - stoich_invariants
    # del_c_iso_red = -del_r_iso_red - del_l_iso_red

    # del_l_iso_red, del_r_iso_red, del_c_iso_red = s_linalg(Svv_sparse_iso_red, S_merge_iso_red, stoichiometric_invariants_iso_red)
else:
    print("Error: numpy vs sympy disparity")
