import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_radm2 = pd.read_csv("../mechanisms/radm2/RADM2_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_radm2, Sp_sparse_radm2, Svv_sparse_radm2 = create_sparse(edge_list_radm2)
Svv_sparse_radm2, init_col_del_radm2 = del_zero_col(Svv_sparse_radm2)
print("computing dimension of nullspace...")
stoichiometric_invariants_radm2 = len(Svv_sparse_radm2.T.nullspace())
Svv_radm2_np = np.array(Svv_sparse_radm2, dtype=float)
rank_Svv_radm2_np = np.linalg.matrix_rank(Svv_radm2_np)
dim_leftnull_radm2 = Svv_sparse_radm2.shape[0] - rank_Svv_radm2_np

print("identifying coproduction columns...")
coproduction_cols_radm2 = create_coproduction(Sr_sparse_radm2)
print("creating symbolic dictionary...")
symbol_dict_radm2 = create_symbols(coproduction_cols_radm2)
print("merging coproduction columns...")
S_merge_radm2, col_del_radm2 = merge_coprod(Sr_sparse_radm2, Sp_sparse_radm2, symbol_dict_radm2, coproduction_cols_radm2, Svv_sparse_radm2)
print("performing linear algebra...")
# rank_radm2 = S_merge_radm2.rank()
# dim_null_radm2 = S_merge_radm2.shape[0] - rank_radm2
# del_l_radm2 = dim_null_radm2 - stoichiometric_invariants_radm2
del_r_radm2 = S_merge_radm2.shape[1] - Svv_sparse_radm2.shape[1]
# del_c_radm2 = -del_r_radm2 - del_l_radm2

rank_list_radm2 = linalg_experiment(S_merge_radm2, num_experiments)
# del_l_radm2 = S_merge_radm2.shape[0] - rank - stoich_invariants
# del_c_radm2 = -del_r_radm2 - del_l_radm2

# del_l_radm2, del_r_radm2, del_c_radm2 = s_linalg(
#     Svv_sparse_radm2,
#     S_merge_radm2,
#     stoichiometric_invariants_radm2
# )
