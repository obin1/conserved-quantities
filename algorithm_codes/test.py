import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment, create_coproduction_2, linalg_experiment_fast, numeric_sparse_matrix_fast_combined

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_iso_full = pd.read_csv("../mechanisms/cri-v2.2/cri_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_iso_full, Sp_sparse_iso_full, Svv_sparse_iso_full = create_sparse(edge_list_iso_full)
print("computing dimension of nullspace...")
# stoichiometric_invariants_iso_full = len(Svv_sparse_iso_full.T.nullspace())
Svv_iso_full_np = np.array(Svv_sparse_iso_full, dtype=float)
rank_Svv_iso_full_np = np.linalg.matrix_rank(Svv_iso_full_np)
dim_leftnull_iso_full = Svv_sparse_iso_full.shape[0] - rank_Svv_iso_full_np

print("identifying coproduction columns...")
coproduction_cols_iso_full = create_coproduction_2(Sr_sparse_iso_full)
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

# rank_list_iso_full = linalg_experiment(S_merge_iso_full, num_experiments)
A = numeric_sparse_matrix_fast_combined(S_merge_iso_full)
rank_list_mcm = linalg_experiment_fast(A, num_experiments, A.shape[0])
# del_l_iso_full = S_merge_iso_full.shape[0] - rank - stoich_invariants
# del_c_iso_full = -del_r_iso_full - del_l_iso_full

# del_l_iso_full, del_r_iso_full, del_c_iso_full = s_linalg(Svv_sparse_iso_full, S_merge_iso_full, stoichiometric_invariants_iso_full)


