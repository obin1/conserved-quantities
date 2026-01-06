import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_form1985 = pd.read_csv("../mechanisms/form1985/form1985_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_form1985, Sp_sparse_form1985, Svv_sparse_form1985 = create_sparse(edge_list_form1985)
print("computing dimension of nullspace...")
stoichiometric_invariants_form1985 = len(Svv_sparse_form1985.T.nullspace())
Svv_form1985_np = np.array(Svv_sparse_form1985, dtype=float)
rank_Svv_form1985_np = np.linalg.matrix_rank(Svv_form1985_np)
dim_leftnull_form1985 = Svv_sparse_form1985.shape[0] - rank_Svv_form1985_np

if stoichiometric_invariants_form1985 == dim_leftnull_form1985:
    print("identifying coproduction columns...")
    coproduction_cols_form1985 = create_coproduction(Sr_sparse_form1985)
    print("creating symbolic dictionary...")
    symbol_dict_form1985 = create_symbols(coproduction_cols_form1985)
    print("merging coproduction columns...")
    S_merge_form1985, col_del_form1985 = merge_coprod(Sr_sparse_form1985, Sp_sparse_form1985, symbol_dict_form1985, coproduction_cols_form1985, Svv_sparse_form1985)
    print("performing linear algebra...")
    # rank_form1985 = S_merge_form1985.rank()
    # dim_null_form1985 = S_merge_form1985.shape[0] - rank_form1985
    # del_l_form1985 = dim_null_form1985 - stoichiometric_invariants_form1985
    del_r_form1985 = S_merge_form1985.shape[1] - Svv_sparse_form1985.shape[1]
    # del_c_form1985 = -del_r_form1985 - del_l_form1985
    
    rank_list_form1985 = linalg_experiment(S_merge_form1985, num_experiments)
    # del_l_form1985 = S_merge_form1985.shape[0] - rank - stoich_invariants
    # del_c_form1985 = -del_r_form1985 - del_l_form1985

    # del_l_form1985, del_r_form1985, del_c_form1985 = s_linalg(Svv_sparse_form19985, S_merge_form1985, stoichiometric_invariants_form1985)
else:
    print("Error: numpy vs sympy disparity")
