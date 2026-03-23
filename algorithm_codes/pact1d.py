import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment
from scipy.linalg import null_space

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_pact1d = pd.read_csv("../mechanisms/pact1d/mech_EdgeList.csv", comment="!")
# this mechanism is from pact-1d halogens, 
# paper https://doi.org/10.1029/2021JD036140
# zenodo https://doi.org/10.5281/zenodo.6045999

print("creating sparse matrices...")
Sr_sparse_pact1d, Sp_sparse_pact1d, Svv_sparse_pact1d = create_sparse(edge_list_pact1d)
Svv_sparse_pact1d, init_col_del_pact1d = del_zero_col(Svv_sparse_pact1d)
print("computing dimension of nullspace...")
Svv_pact1d_np = np.array(Svv_sparse_pact1d, dtype=float)
rank_Svv_pact1d_np = np.linalg.matrix_rank(Svv_pact1d_np)
dim_leftnull_pact1d = Svv_sparse_pact1d.shape[0] - rank_Svv_pact1d_np
pact1d_null_np = null_space(Svv_pact1d_np.T)

print("identifying coproduction columns...")
coproduction_cols_pact1d = create_coproduction(Sr_sparse_pact1d)
print("creating symbolic dictionary...")
symbol_dict_pact1d = create_symbols(coproduction_cols_pact1d)
print("merging coproduction columns...")
S_merge_pact1d, col_del_pact1d = merge_coprod(Sr_sparse_pact1d, Sp_sparse_pact1d, symbol_dict_pact1d, coproduction_cols_pact1d, Svv_sparse_pact1d)
print("performing linear algebra...")
del_r_pact1d = S_merge_pact1d.shape[1] - Svv_sparse_pact1d.shape[1]

rank_list_pact1d = linalg_experiment(S_merge_pact1d, num_experiments)
