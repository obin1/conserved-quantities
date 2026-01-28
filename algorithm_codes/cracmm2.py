import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cracmm = pd.read_csv("../mechanisms/cracmm2/cracmm2_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cracmm, Sp_sparse_cracmm, Svv_sparse_cracmm = create_sparse(edge_list_cracmm)
Svv_sparse_cracmm, init_col_del_cracmm = del_zero_col(Svv_sparse_cracmm)
print("computing dimension of nullspace...")
stoichiometric_invariants_cracmm = len(Svv_sparse_cracmm.T.nullspace())
Svv_cracmm_np = np.array(Svv_sparse_cracmm, dtype=float)
rank_Svv_cracmm_np = np.linalg.matrix_rank(Svv_cracmm_np)
dim_leftnull_cracmm = Svv_sparse_cracmm.shape[0] - rank_Svv_cracmm_np

print("identifying coproduction columns...")
coproduction_cols_cracmm = create_coproduction(Sr_sparse_cracmm)
print("creating symbolic dictionary...")
symbol_dict_cracmm = create_symbols(coproduction_cols_cracmm)
print("merging coproduction columns...")
S_merge_cracmm, col_del_cracmm = merge_coprod(
    Sr_sparse_cracmm,
    Sp_sparse_cracmm,
    symbol_dict_cracmm,
    coproduction_cols_cracmm,
    Svv_sparse_cracmm
)
print("performing linear algebra...")
# rank_cracmm = S_merge_cracmm.rank()
# dim_null_cracmm = S_merge_cracmm.shape[0] - rank_cracmm
# del_l_cracmm = dim_null_cracmm - stoichiometric_invariants_cracmm
del_r_cracmm = S_merge_cracmm.shape[1] - Svv_sparse_cracmm.shape[1]
# del_c_cracmm = -del_r_cracmm - del_l_cracmm
    
rank_list_cracmm = linalg_experiment(S_merge_cracmm, num_experiments)
# del_l_cracmm = S_merge_cracmm.shape[0] - rank - stoich_invariants
# del_c_cracmm = -del_r_cracmm - del_l_cracmm

# del_l_cracmm, del_r_cracmm, del_c_cracmm = s_linalg(
#     Svv_sparse_cracmm,
#     S_merge_cracmm,
#     stoichiometric_invariants_cracmm
# )
