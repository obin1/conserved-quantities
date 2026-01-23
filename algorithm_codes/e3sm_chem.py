import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment
from scipy.linalg import null_space

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_chemuci = pd.read_csv("../mechanisms/e3sm-chem/chemuci_EdgeList.csv", comment="!")
# this mechanism is from pact-1d halogens, 
# paper https://doi.org/10.1029/2021JD036140
# zenodo https://doi.org/10.5281/zenodo.6045999

print("creating sparse matrices...")
Sr_sparse_chemuci, Sp_sparse_chemuci, Svv_sparse_chemuci = create_sparse(edge_list_chemuci)
print("computing dimension of nullspace...")
# chemuci_null = Svv_sparse_chemuci.T.nullspace()
# stoichiometric_invariants_chemuci = len(chemuci_null)

Svv_chemuci_np = np.array(Svv_sparse_chemuci, dtype=float)
rank_Svv_chemuci_np = np.linalg.matrix_rank(Svv_chemuci_np)
dim_leftnull_chemuci = Svv_sparse_chemuci.shape[0] - rank_Svv_chemuci_np
chemuci_null_np = null_space(Svv_chemuci_np.T)

print("identifying coproduction columns...")
coproduction_cols_chemuci = create_coproduction(Sr_sparse_chemuci)
print("creating symbolic dictionary...")
symbol_dict_chemuci = create_symbols(coproduction_cols_chemuci)
print("merging coproduction columns...")
S_merge_chemuci, col_del_chemuci = merge_coprod(Sr_sparse_chemuci, Sp_sparse_chemuci, symbol_dict_chemuci, coproduction_cols_chemuci, Svv_sparse_chemuci)
print("performing linear algebra...")
# rank_chemuci = S_merge_chemuci.rank()
# dim_null_chemuci = S_merge_chemuci.shape[0] - rank_chemuci
# del_l_chemuci = dim_null_chemuci - stoichiometric_invariants_chemuci
del_r_chemuci = S_merge_chemuci.shape[1] - Svv_sparse_chemuci.shape[1]
# del_c_chemuci = -del_r_chemuci - del_l_chemuci

rank_list_chemuci = linalg_experiment(S_merge_chemuci, num_experiments)
# del_l = S_merge_chemuci.shape[0] - rank - stoich_invariants
# del_c_chemuci = -del_r_chemuci - del_l_chemuci