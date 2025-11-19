import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

edge_list_cri = pd.read_csv("../mechanisms/cri-v2.2/cri_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cri, Sp_sparse_cri, Svv_sparse_cri = create_sparse(edge_list_cri)
print("computing dimension of nullspace...")
stoichiometric_invariants_cri = len(Svv_sparse_cri.T.nullspace())
Svv_cri_np = np.array(Svv_sparse_cri, dtype=float)
rank_Svv_cri_np = np.linalg.matrix_rank(Svv_cri_np)
dim_leftnull_cri = Svv_sparse_cri.shape[0] - rank_Svv_cri_np

if stoichiometric_invariants_cri == dim_leftnull_cri:
    print("identifying coproduction columns...")
    coproduction_cols_cri = create_coproduction(Sr_sparse_cri)
    print("creating symbolic dictionary...")
    symbol_dict_cri = create_symbols(coproduction_cols_cri)
    print("merging coproduction columns...")
    S_merge_cri, col_del_cri = merge_coprod(Sr_sparse_cri, Sp_sparse_cri, symbol_dict_cri, coproduction_cols_cri, Svv_sparse_cri)
    print("performing linear algebra...")
    # rank_cri = S_merge_cri.rank()
    # dim_null_cri = S_merge_cri.shape[0] - rank_cri
    # del_l_cri = dim_null_cri - stoichiometric_invariants_cri
    del_r_cri = S_merge_cri.shape[1] - Svv_sparse_cri.shape[1]
    # del_c_cri = -del_r_cri - del_l_cri

    rank_list_cri = linalg_experiment(S_merge_cri)
    # del_l = S_merge_cri.shape[0] - rank - stoich_invariants
    # del_c_cri = -del_r_cri - del_l_cri

    # del_l_cri, del_r_cri, del_c_cri = s_linalg(Svv_sparse_cri, S_merge_cri, stoichiometric_invariants_cri)
else:
    print("Error: numpy vs sympy disparity")

#%%

print("identifying coproduction columns...")
coproduction_cols_cri = create_coproduction(Sr_sparse_cri)
print("creating symbolic dictionary...")
symbol_dict_cri = create_symbols(coproduction_cols_cri)
print("merging coproduction columns...")
S_merge_cri, col_del_cri = merge_coprod(Sr_sparse_cri, Sp_sparse_cri, symbol_dict_cri, coproduction_cols_cri, Svv_sparse_cri)
print("performing linear algebra...")
# rank_cri = S_merge_cri.rank()
# dim_null_cri = S_merge_cri.shape[0] - rank_cri
# del_l_cri = dim_null_cri - stoichiometric_invariants_cri
del_r_cri = S_merge_cri.shape[1] - Svv_sparse_cri.shape[1]
# del_c_cri = -del_r_cri - del_l_cri

rank_list_cri = linalg_experiment(S_merge_cri)
# del_l = S_merge_cri.shape[0] - rank - stoich_invariants
# del_c_cri = -del_r_cri - del_l_cri
# %%
