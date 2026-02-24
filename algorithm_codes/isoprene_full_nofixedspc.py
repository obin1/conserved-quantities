import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, numeric_sparse_matrix_fast_combined, linalg_experiment_fast

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_iso_vo = pd.read_csv("../mechanisms/isoprene_full_v5_nofixedspc/isoprene_full_v5_nofixedspc_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_iso_vo, Sp_sparse_iso_vo, Svv_sparse_iso_vo = create_sparse(edge_list_iso_vo)
Svv_sparse_iso_vo, init_col_del_iso_vo = del_zero_col(Svv_sparse_iso_vo)
print("computing dimension of nullspace...")
stoichiometric_invariants_iso_vo = len(Svv_sparse_iso_vo.T.nullspace())
Svv_iso_vo_np = np.array(Svv_sparse_iso_vo, dtype=float)
rank_Svv_iso_vo_np = np.linalg.matrix_rank(Svv_iso_vo_np)
dim_leftnull_iso_vo = Svv_sparse_iso_vo.shape[0] - rank_Svv_iso_vo_np

print("identifying coproduction columns...")
coproduction_cols_iso_vo = create_coproduction(Sr_sparse_iso_vo)
print("creating symbolic dictionary...")
symbol_dict_iso_vo = create_symbols(coproduction_cols_iso_vo)
print("merging coproduction columns...")
S_merge_iso_vo, col_del_iso_vo = merge_coprod(Sr_sparse_iso_vo, Sp_sparse_iso_vo, symbol_dict_iso_vo, coproduction_cols_iso_vo, Svv_sparse_iso_vo)
print("performing linear algebra...")
# rank_iso_vo = S_merge_iso_vo.rank()
# dim_null_iso_vo = S_merge_iso_vo.shape[0] - rank_iso_vo
# del_l_iso_vo = dim_null_iso_vo - stoichiometric_invariants_iso_vo
del_r_iso_vo = S_merge_iso_vo.shape[1] - Svv_sparse_iso_vo.shape[1]
# del_c_iso_vo = -del_r_iso_vo - del_l_iso_vo

rank_list_iso_vo = linalg_experiment(S_merge_iso_vo, num_experiments)
# del_l_iso_vo = S_merge_iso_vo.shape[0] - rank - stoich_invariants
# del_c_iso_vo = -del_r_iso_vo - del_l_iso_vo

A = numeric_sparse_matrix_fast_combined(S_merge_iso_vo)
rank_list_iso_vo_2 = linalg_experiment_fast(A, num_experiments, A.shape[0])

# del_l_iso_vo, del_r_iso_vo, del_c_iso_vo = s_linalg(Svv_sparse_iso_vo, S_merge_iso_vo, stoichiometric_invariants_iso_vo)

