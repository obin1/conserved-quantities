#%% 
import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

#%%
SEED = 42
np.random.seed(SEED)

num_experiments = 10

pd.set_option('display.max_columns', None)
# edge_list_logan = pd.read_csv("../mechanisms/logan/log81_EdgeList_withloss.csv", comment="!")
# edge_list_logan = pd.read_csv("../mechanisms/logan/log81_EdgeList.csv", comment="!")
edge_list_logan = pd.read_csv("../mechanisms/logan/log81_EdgeList_withN2.csv", comment="!")

print("creating sparse matrices...")
Sr_sparse_logan, Sp_sparse_logan, Svv_sparse_logan = create_sparse(edge_list_logan)
Svv_sparse_logan, init_col_del_logan = del_zero_col(Svv_sparse_logan)
print("computing dimension of nullspace...")
# stoichiometric_invariants_logan = len(Svv_sparse_logan.T.nullspace())
Svv_logan_np = np.array(Svv_sparse_logan, dtype=float)
rank_Svv_logan_np = np.linalg.matrix_rank(Svv_logan_np)
dim_leftnull_logan = Svv_sparse_logan.shape[0] - rank_Svv_logan_np


print("identifying coproduction columns...")
coproduction_cols_logan = create_coproduction(Sr_sparse_logan)
print("creating symbolic dictionary...")
symbol_dict_logan = create_symbols(coproduction_cols_logan)
print("merging coproduction columns...")
S_merge_logan, col_del_logan = merge_coprod(Sr_sparse_logan, Sp_sparse_logan, symbol_dict_logan, coproduction_cols_logan, Svv_sparse_logan)
print("performing linear algebra...")

del_r_logan = S_merge_logan.shape[1] - Svv_sparse_logan.shape[1]
    
rank_list_logan = linalg_experiment(S_merge_logan, num_experiments)
N = S_merge_logan.T.nullspace()

#%%
# for loop over edgelist to get a species index mapping dict. species can either be in from or to columns
# if it begins with "R" then it is a reaction, else it is a species
species_index = {}
for index, row in edge_list_logan.iterrows():
    from_node = row['from']
    to_node = row['to']
    if not from_node.startswith("R"):
        if from_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = from_node
    if not to_node.startswith("R"):
        if to_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = to_node

def get_species_in_null_vector(null_vector, species_index):
    species_involved = []
    for i in range(len(null_vector)):
        if abs(null_vector[i]) != 0: 
            species_name = species_index.get(i + 1, None)  
            if species_name:
                coeff = null_vector[i]
                species_involved.append(f"{coeff}*{species_name}")
    equation = " + ".join(species_involved)
    return equation
first_null_vector = N[1]
species_in_first_null_vector = get_species_in_null_vector(first_null_vector, species_index)
print(species_in_first_null_vector)