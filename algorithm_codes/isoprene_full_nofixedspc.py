import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_iso_full = pd.read_csv("../mechanisms/isoprene_full_v5_nofixedspc/isoprene_full_v5_nofixedspc_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_iso_full, Sp_sparse_iso_full, Svv_sparse_iso_full = create_sparse(edge_list_iso_full)
Svv_sparse_iso_full, init_col_del_iso_full = del_zero_col(Svv_sparse_iso_full)
print("computing dimension of nullspace...")
stoichiometric_invariants_iso_full = len(Svv_sparse_iso_full.T.nullspace())
Svv_iso_full_np = np.array(Svv_sparse_iso_full, dtype=float)
rank_Svv_iso_full_np = np.linalg.matrix_rank(Svv_iso_full_np)
dim_leftnull_iso_full = Svv_sparse_iso_full.shape[0] - rank_Svv_iso_full_np

if stoichiometric_invariants_iso_full == dim_leftnull_iso_full:
    print("identifying coproduction columns...")
    coproduction_cols_iso_full = create_coproduction(Sr_sparse_iso_full)
    print("creating symbolic dictionary...")
    symbol_dict_iso_full = create_symbols(coproduction_cols_iso_full)
    print("merging coproduction columns...")
    S_merge_iso_full, col_del_iso_full = merge_coprod(Sr_sparse_iso_full, Sp_sparse_iso_full, symbol_dict_iso_full, coproduction_cols_iso_full, Svv_sparse_iso_full)
    print("performing linear algebra...")
    # rank_iso_full = S_merge_iso_full.rank()
    # dim_null_iso_full = S_merge_iso_full.shape[0] - rank_iso_full
    # del_l_iso_full = dim_null_iso_full - stoichiometric_invariants_iso_full
    del_r_iso_full = S_merge_iso_full.shape[1] - Svv_sparse_iso_full.shape[1]
    # del_c_iso_full = -del_r_iso_full - del_l_iso_full
    
    rank_list_iso_full = linalg_experiment(S_merge_iso_full, num_experiments)
    # del_l_iso_full = S_merge_iso_full.shape[0] - rank - stoich_invariants
    # del_c_iso_full = -del_r_iso_full - del_l_iso_full

    # del_l_iso_full, del_r_iso_full, del_c_iso_full = s_linalg(Svv_sparse_iso_full, S_merge_iso_full, stoichiometric_invariants_iso_full)
else:
    print("Error: numpy vs sympy disparity")
