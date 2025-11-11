#%%

import pandas as pd
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

#%%

edge_list_gc = pd.read_csv("../mechanisms/geos-chem-v14/gckpp_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_gc, Sp_sparse_gc, Svv_sparse_gc = create_sparse(edge_list_gc)
print("computing dimension of nullspace...")
stoichiometric_invariants_gc = len(Svv_sparse_gc.T.nullspace())
print("identifying coproduction columns...")
coproduction_cols_gc = create_coproduction(Sr_sparse_gc)
print("creating symbolic dictionary...")
symbol_dict_gc = create_symbols(coproduction_cols_gc)
print("merging coproduction columns...")
S_merge_gc = merge_coprod(Sr_sparse_gc, Sp_sparse_gc, symbol_dict_gc, coproduction_cols_gc, Svv_sparse_gc)
print("performing linear algebra...")
rank_gc = S_merge_gc.rank()
dim_null_gc = S_merge_gc.shape[0] - rank_gc
del_l_gc = dim_null_gc - stoichiometric_invariants_gc
del_r_gc = S_merge_gc.shape[1] - Svv_sparse_gc.shape[1]
del_c_gc = -del_r_gc - del_l_gc
# del_l_gc, del_r_gc, del_c_gc = s_linalg(Svv_sparse_gc, S_merge_gc, stoichiometric_invariants_gc)

