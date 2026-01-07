import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cbmz = pd.read_csv("../mechanisms/cbmz/CBMZ_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cbmz, Sp_sparse_cbmz, Svv_sparse_cbmz = create_sparse(edge_list_cbmz)
print("computing dimension of nullspace...")
stoichiometric_invariants_cbmz = len(Svv_sparse_cbmz.T.nullspace())
Svv_cbmz_np = np.array(Svv_sparse_cbmz, dtype=float)
rank_Svv_cbmz_np = np.linalg.matrix_rank(Svv_cbmz_np)
dim_leftnull_cbmz = Svv_sparse_cbmz.shape[0] - rank_Svv_cbmz_np

if stoichiometric_invariants_cbmz == dim_leftnull_cbmz:
    print("identifying coproduction columns...")
    coproduction_cols_cbmz = create_coproduction(Sr_sparse_cbmz)
    print("creating symbolic dictionary...")
    symbol_dict_cbmz = create_symbols(coproduction_cols_cbmz)
    print("merging coproduction columns...")
    S_merge_cbmz, col_del_cbmz = merge_coprod(Sr_sparse_cbmz, Sp_sparse_cbmz, symbol_dict_cbmz, coproduction_cols_cbmz, Svv_sparse_cbmz)
    print("performing linear algebra...")
    # rank_cbmz = S_merge_cbmz.rank()
    # dim_null_cbmz = S_merge_cbmz.shape[0] - rank_cbmz
    # del_l_cbmz = dim_null_cbmz - stoichiometric_invariants_cbmz
    del_r_cbmz = S_merge_cbmz.shape[1] - Svv_sparse_cbmz.shape[1]
    # del_c_cbmz = -del_r_cbmz - del_l_cbmz
    
    rank_list_cbmz = linalg_experiment(S_merge_cbmz, num_experiments)
    # del_l_cbmz = S_merge_cbmz.shape[0] - rank - stoich_invariants
    # del_c_cbmz = -del_r_cbmz - del_l_cbmz

    # del_l_cbmz, del_r_cbmz, del_c_cbmz = s_linalg(
    #     Svv_sparse_cbmz,
    #     S_merge_cbmz,
    #     stoichiometric_invariants_cbmz
    # )
else:
    print("Error: numpy vs sympy disparity")