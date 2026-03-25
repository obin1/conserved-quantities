import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, del_zero_col, create_coproduction, create_symbols, merge_coprod, s_linalg, linalg_experiment

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Set number of experiments for linear algebra step
num_experiments = 10

# Define stoichiometric matrix
Svv_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [-2, -1,  0,  0, -1,  1, -1,  1,  1,  0], # O2
    [ 1, -2, -2,  0,  0, -1,  0, -2,  0,  0], # PhCH2O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 1,  2,  0,  0,  1, -1,  1,  0, -2,  0], # HO2
    [ 0,  2,  1, -1,  1,  0,  1,  0,  0,  0], # PhCHO
    [ 0,  0,  1,  0, -1,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0, -1,  0,  1,  0,  0,  2], # OH
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  1, -1,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1], # H2O2
])
# Convert to sparse matrix
Svv_sparse_d = sp.SparseMatrix(Svv_d)

# Define reactant matrix
Sr_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [-2, -1,  0,  0, -1,  0, -1,  0,  0,  0], # O2
    [ 0, -2, -2,  0,  0, -1,  0, -2,  0,  0], # PhCH2O2
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 0,  0,  0,  0,  0, -1,  0,  0, -2,  0], # HO2
    [ 0,  0,  0, -1,  0,  0,  0,  0,  0,  0], # PhCHO
    [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # OH
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # H2O2
])
# Convert to sparse matrix
Sr_sparse_d = sp.SparseMatrix(Sr_d)

# Define product matrix
Sp_d = sp.Matrix([
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [ 0,  0,  0,  0,  0,  1,  0,  1,  1,  0], # O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2O2
    [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # CO
    [ 1,  2,  0,  0,  1,  0,  1,  0,  0,  0], # HO2
    [ 0,  2,  1,  0,  1,  0,  1,  0,  0,  0], # PhCHO
    [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  0], # PhCH2OH
    [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  2], # OH
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # H2O
    [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0], # PhCH2OOH
    [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # PhCH2O2CH2Ph
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # H2O2
])
# Convert to sparse matrix
Sp_sparse_d = sp.SparseMatrix(Sp_d)

# Compute dimension of left null space in two ways 
# (SymPy nullspace operations sometimes result in a different number due to precision errors)
print("computing dimension of nullspace...")
# Compute dimension of left null space using SymPy nullspace operation
stoichiometric_invariants_d = len(Svv_sparse_d.T.nullspace())
# Convert stoichiometric SymPy matrix to NumPy matrix
Svv_d_np = np.array(Svv_sparse_d, dtype=float)
# Compute rank of NumPy matrix
rank_Svv_d_np = np.linalg.matrix_rank(Svv_d_np)
# Compute dimension of left null space by subtracting rank from the number of rows
dim_leftnull_d = Svv_sparse_d.shape[0] - rank_Svv_d_np

# Identify sets of coproducing reactions
print("identifying coproduction columns...")
coproduction_cols_d = create_coproduction(Sr_sparse_d)
# Create the dictionary of symbols for coproducing groups
print("creating symbolic dictionary...")
symbol_dict_d = create_symbols(coproduction_cols_d)
# Perform the column "merging" operation on branching reaction groups
print("merging coproduction columns...")
S_merge_d, col_del_d = merge_coprod(Sr_sparse_d, Sp_sparse_d, symbol_dict_d, coproduction_cols_d, Svv_sparse_d)
print("performing linear algebra...")
# Compute the number of reactions lost due to merging
del_r_d = S_merge_d.shape[1] - Svv_sparse_d.shape[1]
# Perform the rank calculation experiment on S_merge
# The number of kinetic invariants = S_merge_d.shape[0] - rank(S_merge) - # stoichiometric invariants
rank_list_d = linalg_experiment(S_merge_d, num_experiments)