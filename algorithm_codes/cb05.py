import pandas as pd
import numpy as np
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

SEED = 42
np.random.seed(SEED)

num_experiments = 10

edge_list_cb05 = pd.read_csv("../mechanisms/cb05/CB05TUCl_EPA_EdgeList.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_cb05, Sp_sparse_cb05, Svv_sparse_cb05 = create_sparse(edge_list_cb05)
print("computing dimension of nullspace...")
stoichiometric_invariants_cb05 = len(Svv_sparse_cb05.T.nullspace())
Svv_cb05_np = np.array(Svv_sparse_cb05, dtype=float)
rank_Svv_cb05_np = np.linalg.matrix_rank(Svv_cb05_np)
dim_leftnull_cb05 = Svv_sparse_cb05.shape[0] - rank_Svv_cb05_np

if stoichiometric_invariants_cb05 == dim_leftnull_cb05:
    print("identifying coproduction columns...")
    coproduction_cols_cb05 = create_coproduction(Sr_sparse_cb05)
    print("creating symbolic dictionary...")
    symbol_dict_cb05 = create_symbols(coproduction_cols_cb05)
    print("merging coproduction columns...")
    S_merge_cb05, col_del_cb05 = merge_coprod(
        Sr_sparse_cb05,
        Sp_sparse_cb05,
        symbol_dict_cb05,
        coproduction_cols_cb05,
        Svv_sparse_cb05
    )
    print("performing linear algebra...")
    # rank_cb05 = S_merge_cb05.rank()
    # dim_null_cb05 = S_merge_cb05.shape[0] - rank_cb05
    # del_l_cb05 = dim_null_cb05 - stoichiometric_invariants_cb05
    del_r_cb05 = S_merge_cb05.shape[1] - Svv_sparse_cb05.shape[1]
    # del_c_cb05 = -del_r_cb05 - del_l_cb05
    
    rank_list_cb05 = linalg_experiment(S_merge_cb05, num_experiments)
    # del_l_cb05 = S_merge_cb05.shape[0] - rank - stoich_invariants
    # del_c_cb05 = -del_r_cb05 - del_l_cb05

    # del_l_cb05, del_r_cb05, del_c_cb05 = s_linalg(
    #     Svv_sparse_cb05,
    #     S_merge_cb05,
    #     stoichiometric_invariants_cb05
    # )
else:
    print("Error: numpy vs sympy disparity")
