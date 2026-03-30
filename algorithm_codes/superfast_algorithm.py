import pandas as pd
import numpy as np
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Read mechanism EdgeList (CSV) into a Pandas DataFrame
edge_list_superfast = pd.read_csv("../mechanisms/superfast/superfast_EdgeList.csv", comment="!")

# Create the sparse stoichiometric, reactant, and product matrices
print("creating sparse matrices...")
Sr_sparse_superfast, Sp_sparse_superfast, Svv_sparse_superfast = create_sparse(edge_list_superfast)
# Remove any 0 columns unrelated to merging
Svv_sparse_superfast, Sp_sparse_superfast, Sr_sparse_superfast, init_col_del_superfast = del_zero_col(Svv_sparse_superfast, Sp_sparse_superfast, Sr_sparse_superfast)
print("computing dimension of nullspace...")
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_superfast_np = np.array(Svv_sparse_superfast, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_superfast_np = np.linalg.matrix_rank(Svv_superfast_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_superfast = Svv_sparse_superfast.shape[0] - rank_Svv_superfast_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_superfast = create_coproduction(Sr_sparse_superfast)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_superfast = create_symbols(coproduction_cols_superfast)
# Perform the column "merging" operation on coproducting groups
print("merging coproduction columns...")
S_merge_superfast, col_del_supefast = merge_coprod(Sr_sparse_superfast, Sp_sparse_superfast, symbol_dict_superfast, coproduction_cols_superfast, Svv_sparse_superfast)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_superfast = S_merge_superfast.shape[1] - Svv_sparse_superfast.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge.shape[0] - rank(S_merge) - # stoichiometric invariants 
rank_list_superfast = linalg_experiment(S_merge_superfast, num_experiments)

# Use SymPy null space operation to identify the single kinetic invariant
N = S_merge_superfast.T.nullspace()

# Loop over the edge list to create a species index mapping dictionary 
# Species can either be in from or to columns
# If it begins with "R" then it is a reaction number, else it is a species
species_index = {}
for index, row in edge_list_superfast.iterrows():
    from_node = row['from']
    to_node = row['to']
    if not from_node.startswith("R"):
        if from_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = from_node
    if not to_node.startswith("R"):
        if to_node not in species_index:
            species_index[row["# species_index (starts from 1)"]] = to_node

# Construct the kinetic invariant vector as a string combination of the involved species, with their respective coefficients
def get_species_in_null_vector(null_vector, species_index):
    species_involved = []
    for i in range(len(null_vector)):
        if abs(null_vector[i]) != 0: 
            species_name = species_index.get(i + 1, None)  
            if species_name:
                coeff = null_vector[i]
                species_involved.append(f"({coeff})*{species_name}")
    equation = " + ".join(species_involved)
    return equation

# Obtain the vector in the left null space of S_merge_logan
first_null_vector = N[0]
species_in_first_null_vector = get_species_in_null_vector(first_null_vector, species_index)
print(species_in_first_null_vector)