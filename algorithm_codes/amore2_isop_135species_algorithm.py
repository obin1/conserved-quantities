import pandas as pd
import numpy as np
from sympy import S
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_amore = pd.read_csv("../mechanisms/amore2_isop_135species/amore2_isop_135species_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_amore, Sp_sparse_amore, Svv_sparse_amore = create_sparse(edge_list_amore)
Svv_sparse_amore, init_col_del_amore = del_zero_col(Svv_sparse_amore)

zero_rows_indices = []
for i in range(Svv_sparse_amore.rows):
    row = Svv_sparse_amore.row(i)
    if all(element == S.Zero for element in row):
        zero_rows_indices.append(i)
for i in reversed(zero_rows_indices):
    Svv_sparse_amore.row_del(i)
    Sr_sparse_amore.row_del(i)
    Sp_sparse_amore.row_del(i)

print("computing dimension of nullspace...")
stoichiometric_invariants_amore = len(Svv_sparse_amore.T.nullspace())
Svv_amore_np = np.array(Svv_sparse_amore, dtype=float)
rank_Svv_amore_np = np.linalg.matrix_rank(Svv_amore_np)
dim_leftnull_amore = Svv_sparse_amore.shape[0] - rank_Svv_amore_np

if stoichiometric_invariants_amore == dim_leftnull_amore:
    print("identifying coproduction columns...")
    coproduction_cols_amore = create_coproduction(Sr_sparse_amore)
    print("creating symbolic dictionary...")
    symbol_dict_amore = create_symbols(coproduction_cols_amore)
    print("merging coproduction columns...")
    S_merge_amore, col_del_amore = merge_coprod(Sr_sparse_amore, Sp_sparse_amore, symbol_dict_amore, coproduction_cols_amore, Svv_sparse_amore)
    print("performing linear algebra...")
    # rank_amore = S_merge_amore.rank()
    # dim_null_amore = S_merge_amore.shape[0] - rank_amore
    # del_l_amore = dim_null_amore - stoichiometric_invariants_amore
    del_r_amore = S_merge_amore.shape[1] - Svv_sparse_amore.shape[1]
    # del_c_amore = -del_r_amore - del_l_amore

    rank_list_amore = linalg_experiment(S_merge_amore, num_experiments)
    # del_l = S_merge_amore.shape[0] - rank - stoich_invariants
    # del_c_amore = -del_r_amore - del_l_amore

    # del_l_amore, del_r_amore, del_c_amore = s_linalg(Svv_sparse_amore, S_merge_amore, stoichiometric_invariants_amore)
else:
    print("Error: numpy vs sympy disparity")

# %%
