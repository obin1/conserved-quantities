import numpy as np
import sympy as sp
# from sympy.matrices.normalforms import DomainMatrix
# from sympy.polys.domains import QQ

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

# ------------------------------------
# Create coproduction function for MCM 
# ------------------------------------
from collections import defaultdict

def create_coproduction_2(Sr):
    n_rows, n_cols = Sr.shape

    groups = defaultdict(list)

    for j in range(n_cols):
        # Get nonzero structure of column j
        col_entries = Sr[:, j].todok()

        # Canonical, hashable signature
        signature = tuple(sorted(col_entries.items()))
        # each item is (row_index, symbolic_value)

        groups[signature].append(j)

    coproduction_cols = [0] * n_cols

    for group in groups.values():
        if len(group) > 1:
            rep = min(group) + 1
            for j in group:
                coproduction_cols[j] = rep

    return coproduction_cols
# ------------------------------------


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
        col_del = []
        for j in reversed(range(S_merge.shape[1])):
            if all(S_merge[i, j] == 0 for i in range(S_merge.shape[0])):
                S_merge.col_del(j)
                col_del.append(j)

        S_merge = S_merge.applyfunc(sp.simplify)

    else:
        S_merge = Svv_sparse
        col_del = []

    return S_merge, col_del


def s_linalg(Svv_sparse, S_merge, stoich_invariants):
    l_vv = stoich_invariants
    l_merge = len(S_merge.T.nullspace())
    del_l = l_merge - l_vv

    # del_c = number of broken cycles
    # -del_r - del_l = del_c
    # reactions in Svv_sparse > reactions in S_merge, so del_r < 0
    del_r = S_merge.shape[1] - Svv_sparse.shape[1]
    del_c = -del_r - del_l

    return del_l, del_r, del_c


def linalg_experiment(S_merge, num_experiments):
    rank_list = []
    for i in range(num_experiments):
        subs_dict = {}
        for key in S_merge.free_symbols:
            subs_dict[key] = np.random.uniform(1, 10)
        subs_matrix = S_merge.subs(subs_dict).evalf()
        subs_matrix_np = np.array(subs_matrix.tolist(), dtype=float)
        rank_list.append(np.linalg.matrix_rank(subs_matrix_np))
    return rank_list


# ---------------------------------
# Linear algebra function for MCM 
# ---------------------------------
def linalg_experiment_fast(A, num_experiments, proj_dim):
    n_rows, n_cols = A.shape
    rank_list = []

    for _ in range(num_experiments):
    
        # 2. Random projection
        R = np.random.randn(n_cols, proj_dim)
        AR = A @ R   # shape: (n_rows, proj_dim)

        # 3. Rank of small matrix
        rank = np.linalg.matrix_rank(AR, tol=1e-8)
        rank_list.append(rank)

    return rank_list
# ---------------------------------


# --------------------------------------------------------------
# Create sparse matrix for faster rank solving function for MCM 
# --------------------------------------------------------------
import scipy
import sympy as sp
import numpy as np
from scipy.sparse import csc_matrix

def numeric_sparse_matrix_fast_combined(S, seed=None):
    """
    Convert a large symbolic SymPy sparse matrix S into a numeric SciPy CSC sparse matrix quickly.
    Combines:
      1. Precomputed DOK entries
      2. Vectorized substitution via lambdify
      3. Bulk float conversion
    """
    if seed is not None:
        np.random.seed(seed)

    print("Precompute DOK entries")
    dok = S.todok()
    keys = list(dok.keys())         # list of (i,j)
    values = list(dok.values())     # list of sympy expressions

    if not values:
        return csc_matrix(S.shape)

    print("Get all free symbols in order")
    free_syms = sorted(list(S.free_symbols), key=lambda s: s.name)

    print("Create fast numeric function using lambdify")
    f = sp.lambdify(free_syms, values, modules='numpy')

    print("Sample random values for symbols")
    vals = np.random.uniform(1.0, 10.0, size=len(free_syms))

    print("Evaluate all nonzero entries at once (vectorized)")
    data = np.array(f(*vals), dtype=np.float64)

    print("Extract row and column indices")
    rows, cols = zip(*keys)

    print("Build CSC sparse matrix")
    return csc_matrix((data, (rows, cols)), shape=S.shape)
# ------------------------------------


# to check if all lost reactions result from coproduction merging 
# (not weird reactions like rxn 61 in SAPRC99)
def check_merge(coproduction_cols):
    x = [i for i in coproduction_cols if i!= 0]
    y = set([i for i in coproduction_cols if i!= 0])
    del_r = -(len(x) - len(y))
    return del_r

# to check EdgeList vs .eqn vs sparse matrix
# np.where(np.array(Svv_sparse[:, n]) != 0)

# check if nullspace of Svv still in nullspace of Smerge
# for v in Svv_sparse.T.nullspace():
    # w = (subs_matrix.T * v).applyfunc(sp.simplify)
    # if not w.is_zero_matrix:
    #     print("False")
    # else: 
    #     print("True")


# a) len(col_del)
# b) len(set(coproduction_cols))
# c) len([x for x in coproduction_cols if x > 0])
# c) - b) + 1 = a)
