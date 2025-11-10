import numpy as np
import sympy as sp
from sympy.matrices.normalforms import DomainMatrix
from sympy.polys.domains import QQ

def create_sparse(edge_list):
    consumed = edge_list[edge_list[" directed stoichiometric value"]<0]
    produced = edge_list[edge_list[" directed stoichiometric value"]>0]
    num_species = edge_list["# species_index (starts from 1)"].max()
    num_reactions = edge_list["reaction_index (starts from 1)"].max()

    Sr = sp.MutableSparseMatrix(num_species, num_reactions, {})
    Sp = sp.MutableSparseMatrix(num_species, num_reactions, {})

    for _, row in consumed.iterrows():
        species_idx = int(row["# species_index (starts from 1)"]) - 1
        reaction_idx = int(row["reaction_index (starts from 1)"]) - 1
        value = float(row[" directed stoichiometric value"])
        Sr[species_idx, reaction_idx] = value

    for _, row in produced.iterrows():
        species_idx = int(row["# species_index (starts from 1)"]) - 1
        reaction_idx = int(row["reaction_index (starts from 1)"]) - 1
        value = float(row[" directed stoichiometric value"])
        Sp[species_idx, reaction_idx] = value

    St = Sr + Sp
    Sr_sparse = sp.SparseMatrix(Sr)
    Sp_sparse = sp.SparseMatrix(Sp)
    Svv_sparse = sp.SparseMatrix(St)

    return Sr_sparse, Sp_sparse, Svv_sparse


def create_coproduction(Sr):
    n_cols = Sr.shape[1]
    coproduction_cols = list(range(1, n_cols + 1))

    for i in range(n_cols):
        col_i = Sr[:, i]
        for j in range(i+1, n_cols):
            col_j = Sr[:, j]
            if col_i == col_j:
                rep = min(coproduction_cols[i], coproduction_cols[j])
                coproduction_cols[j] = rep
    
    for idx in range(n_cols):
        if coproduction_cols.count(coproduction_cols[idx]) == 1:
            coproduction_cols[idx] = 0

    return coproduction_cols


def create_symbols(coproduction_cols):
    if coproduction_cols != []:
        symbol_dict = {}
        coprod_array = np.array(coproduction_cols)
        for num in set(coproduction_cols):
            if num != 0:
                indices = np.where(coprod_array == num)[0] + 1
                symbols = [sp.symbols(f"k{i}") for i in indices]
                total = sum(symbols)
                for item, symbol in zip(indices, symbols):
                    symbol_dict[item.item()] = symbol / total
    else:
        symbol_dict = []

    return symbol_dict


def merge_coprod(Sr, Sp, symbol_dict, coproduction_cols, Svv_sparse):
    if coproduction_cols != []:
        for key, val in symbol_dict.items():
            Sp[:, key-1] = Sp[:, key-1]*val

        # initialize empty matrix 
        n_rows, n_cols = Sr.shape
        S_merge = sp.SparseMatrix.zeros(n_rows, n_cols)

        # fill in the non-coproduction columns
        for idx, label in enumerate(coproduction_cols):
            if label == 0:
                S_merge[:, idx] = Sr[:, idx] + Sp[:, idx]

        # for group in coproduction_cols, perform the merge, take note of first index in the group
        # place merged column in the first index position
        coprod_array = np.array(coproduction_cols)
        unique_groups = set(c for c in coproduction_cols if c != 0)

        for num in unique_groups:
            indices = np.where(coprod_array == num)[0] 
            cols_list = [Sp[:, idx] for idx in indices]
            sum_col = sum(cols_list, sp.SparseMatrix.zeros(n_rows, 1) + Sr[:, indices[0]])
            S_merge[:, indices[0]] = sum_col
        
        # drop columns with only 0, avoid shifting of indices by reversing range
        for j in reversed(range(S_merge.shape[1])):
            if all(S_merge[i, j] == 0 for i in range(S_merge.shape[0])):
                S_merge.col_del(j)

        S_merge = S_merge.applyfunc(sp.simplify)

    else:
        S_merge = Svv_sparse

    return S_merge


def s_linalg(Svv_sparse, S_merge, stoich_invariants):
    l_vv = stoich_invariants
    print("getting dim of nullspace of S_merge...")
    l_merge = len(S_merge.T.nullspace())
    del_l = l_merge - l_vv

    # del_c = number of broken cycles
    # -del_r - del_l = del_c
    # reactions in Svv_sparse > reactions in S_merge, so del_r < 0
    del_r = S_merge.shape[1] - Svv_sparse.shape[1]
    del_c = -del_r - del_l

    return del_l, del_r, del_c


