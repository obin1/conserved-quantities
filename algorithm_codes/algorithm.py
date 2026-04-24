import numpy as np
import sympy as sp
from collections import defaultdict

"""
Creates reactant, product, and full stoichiometric sparse matrices using mechanism edge lists.
Args:
    edge_list (Pandas DataFrame read from csv): Edge list of species-reaction matrix

Returns:
    Sr_sparse (sparse SymPy matrix:) Consumption sparse matrix
    Sp_sparse (sparse SymPy matrix): Production sparse matrix
    Svv_sparse (sparse SymPy matrix): Stoichiometric sparse matrix (net stoichiometry)
"""
def create_sparse(edge_list):
    # We separate species consumption and species production.
    consumed = edge_list[edge_list[" directed stoichiometric value"]<0]
    produced = edge_list[edge_list[" directed stoichiometric value"]>0]

    # Identify total number of species and reactions in the mechanism.
    num_species = edge_list["# species_index (starts from 1)"].max()
    num_reactions = edge_list["reaction_index (starts from 1)"].max()

    # Initialize 2 empty sparse SymPy matrices: one for consumption (Sr) and one for production (Sp).
    Sr = sp.MutableSparseMatrix(num_species, num_reactions, {})
    Sp = sp.MutableSparseMatrix(num_species, num_reactions, {})

    # Populate Sr with data from the edge list.
    for _, row in consumed.iterrows():
        species_idx = int(row["# species_index (starts from 1)"]) - 1
        reaction_idx = int(row["reaction_index (starts from 1)"]) - 1
        value = float(row[" directed stoichiometric value"])
        Sr[species_idx, reaction_idx] = value

    # Populate Sp with data from the edge list.
    for _, row in produced.iterrows():
        species_idx = int(row["# species_index (starts from 1)"]) - 1
        reaction_idx = int(row["reaction_index (starts from 1)"]) - 1
        value = float(row[" directed stoichiometric value"])
        Sp[species_idx, reaction_idx] = value

    # Combine St and Sp to obtain the full stoichiometric matrix.
    St = Sr + Sp

    Sr_sparse = sp.SparseMatrix(Sr)
    Sp_sparse = sp.SparseMatrix(Sp)
    Svv_sparse = sp.SparseMatrix(St)

    return Sr_sparse, Sp_sparse, Svv_sparse

"""
Some reactions both consume and produce the same species, resulting in 0 net change.
We delete these empty reaction columns from the stoichiometric matrix.

Args:
    Svv (SymPy sparse matrix): full stoichiometric matrix of a mechanism (net stoichiometry)

Returns:
    Svv (SymPy sparse matrix): full stoichiometric matrix with zero columns deleted (net stoichiometry)
    col_del (list): list of the deleted columns 
"""
def del_zero_col(Svv, Sp, Sr):
    # Initialize an empty list
    col_del = []

    # Loop through the columns of Svv
    for j in reversed(range(Svv.shape[1])):
        # Check if column is all 0, then append to col_del if yes
        if all(Svv[i, j] == 0 for i in range(Svv.shape[0])):
            Svv.col_del(j)
            Sp.col_del(j)
            Sr.col_del(j)
            col_del.append(j)

    return Svv, Sp, Sr, col_del

"""
Identify which reactions participate in coproduction and group them together by numerical identifiers.

Args:
    Sr (SymPy sparse matrix): reactant matrix of a mechanism

Returns:
    coproduction_cols: list of column indices representing coproduction groupings 
                        (e.g. 1 indicates all reactions that coproduce with reaction 1; 
                              0 indicates that a reaction does not coproduce with any other) 
"""
def create_coproduction(Sr):
    # Get the number of rows and columns in the matrix
    n_rows, n_cols = Sr.shape

    # Initialize an empty defaultdict expecting lists as values
    groups = defaultdict(list)

    for j in range(n_cols):
        # Extract column j as a dictionary of nonzero entries
        # Format: {(row_index, col_index): value}
        col_entries = Sr[:, j].todok()
        # Convert dictionary items into a sorted tuple
        signature = tuple(sorted(col_entries.items()))
        # In the dictionary, group columns that have the same signature
        groups[signature].append(j)

    # Initialize output list with 0's
    coproduction_cols = [0] * n_cols

    for group in groups.values():
        # Only care about groups with more than one column (duplicates)
        if len(group) > 1:
            # Choose a representative column index (smallest index + 1 for 1-based indexing)
            rep = min(group) + 1
            # Assign this representative to all columns in the group
            for j in group:
                coproduction_cols[j] = rep

    return coproduction_cols

"""
Create symbolic variables for each reaction in a coproduction group

Args:
    coproduction_cols: list of column indices representing coproduction groupings 

Returns:
    symbol_dict: dictionary of fractional kinetic rates (symbolic variables for SymPy matrices)
                (e.g. if reactions 1 and 2 coproduct, the dictionary will contain 
                 {1: k1/(k1 + k2), 2: k2/(k1 + k2)}) 
"""
def create_symbols(coproduction_cols):
    if coproduction_cols != []:
        # Initialize an empty dictionary
        symbol_dict = {}
        # Turn coproduction_cols (list) into a NumPy array
        coprod_array = np.array(coproduction_cols)
        # Obtain only unique indices in coproduction_cols
        for num in set(coproduction_cols):
            # Exclude index 0 since that indicates no coproduction
            if num != 0:
                # Identify the array positions of reactions that share coproduction grouping index (1-based counting)
                indices = np.where(coprod_array == num)[0] + 1
                # Create symbolic variables for each participating reaction
                symbols = [sp.symbols(f"k{i}") for i in indices]
                # Create a symbolic sum to be the denominator of the fractional kinetic rates
                total = sum(symbols)
                # Create symbolic fraction and add it to the dictionary, with reaction index as key
                for item, symbol in zip(indices, symbols):
                    symbol_dict[item.item()] = symbol / total
    else:
        symbol_dict = []

    return symbol_dict


"""
Merge coproducing reactions in the stoichiometric matrix

Args:
    Sr: (SymPy sparse matrix): reactant matrix of a mechanism
    Sp: (SymPy sparse matrix): product matrix of a mechanism
    symbol_dict: dictionary of symbolic fractional kinetic rate variables for each reaction in a coproduction group
    coproduction_cols: list of column indices representing coproduction groupings 
    Svv_sparse (SymPy sparse matrix): full stoichiometric matrix of a mechanism (net stoichiometry)

Returns:
    S_merge (SymPy matrix): stoichiometric matrix after merging coproducing reactions
    col_del: list of indices of reactions that were deleted due to merging
"""
def merge_coprod(Sr, Sp, symbol_dict, coproduction_cols, Svv_sparse):
    if coproduction_cols != []:
        # Loop over reaction index: symbolic fractional kinetic rate pairings in symbol_dict
        for key, val in symbol_dict.items():
            # Multiply product coefficients by new fractional kinetic rates for each reaction
            Sp[:, key-1] = Sp[:, key-1]*val

        # Initialize an empty matrix 
        n_rows, n_cols = Sr.shape
        S_merge = sp.SparseMatrix.zeros(n_rows, n_cols)

        # If a reaction does not coproduct, populate that column in S_merge using Sr and Sp values for that index
        for idx, label in enumerate(coproduction_cols):
            if label == 0:
                S_merge[:, idx] = Sr[:, idx] + Sp[:, idx]

        coprod_array = np.array(coproduction_cols)
        unique_groups = set(c for c in coproduction_cols if c != 0)

        # Loop over unique coproduction reaction indices
        for num in unique_groups:
            # Locate other reactions that share the coproduction index
            indices = np.where(coprod_array == num)[0] 
            # Group and sum all participating reaction columns in Sp (with new fractional rates) for this coproduction index
            cols_list = [Sp[:, idx] for idx in indices]
            # Add the corresponding reaction columns in Sr
            sum_col = sum(cols_list, sp.SparseMatrix.zeros(n_rows, 1) + Sr[:, indices[0]])
            # Place new merged column in S_merge at the smallest reaction index in that coproduction group
            S_merge[:, indices[0]] = sum_col
        
        # If columns only contain zeroes, they were merged into a single column
        # Drop these columns
        col_del = []
        # Avoid shifting of indices when deleting columns by looping over the columns in reverse order
        for j in reversed(range(S_merge.shape[1])):
            if all(S_merge[i, j] == 0 for i in range(S_merge.shape[0])):
                S_merge.col_del(j)
                col_del.append(j)

        # Simplify any arithmetic caused by multiplication and addition of symbolic variables
        S_merge = S_merge.applyfunc(sp.simplify)

    # If there was no coproduction in the mechanism, S_merge is equal to Svv
    else:
        S_merge = Svv_sparse
        col_del = []

    return S_merge, col_del

"""
Perform linear algebra to identify # kinetic invariants, # lost reactions, # broken null cycles

Args:
    Svv_sparse (SymPy sparse matrix): full stoichiometric matrix of a mechanism (net stoichiometry)
    S_merge (SymPy matrix): stoichiometric matrix after merging coproducing reactions
    stoich_invariants (int): the number of stoichiometric (non-kinetic) invariants of the original system

Returns:
    del_l: number of kinetic invariants
    del_r: number of reactions lost from merging
    del_c: number of broken null cycles
"""
def s_linalg(Svv_sparse, S_merge, stoich_invariants):
    l_vv = stoich_invariants
    # Obtain the dimensionality of the left null space of S_merge
    l_merge = len(S_merge.T.nullspace())
    # The number of kinetic reactions = dim(left null space of S_merge) - # stoichiometric invariants
    del_l = l_merge - l_vv

    # Obtain the number of columns (reactions) lost between Svv_sparse and S_merge due to merging
    del_r = S_merge.shape[1] - Svv_sparse.shape[1]
    # The number of lost reaction minus the number of kinetic invariants equals the number of broken null cycles
    del_c = -del_r - del_l

    return del_l, del_r, del_c

"""
For larger mechanisms, instead of performing linear algebra on symbolic matrices, substitute uniformly random 
values into the symbols and evaluate the rank of the matrix with substitutions. Repeat for desired number 
of experiments.

Args:
    S_merge (SymPy matrix): stoichiometric matrix after merging coproducing reactions
    num_experiments (int): the number of desired experiments

Returns:
    rank_list (list): a list containing the ranks of the matrix with substituted values, repeated per experiment
"""
def linalg_experiment(S_merge, num_experiments):
    # Initialize an empty list
    rank_list = []
    # Loop over the number of desired experiments.
    for i in range(num_experiments):
        # Initialize an empty dictionary
        subs_dict = {}
        # Loop over the SymPy symbols present in the matrix
        for key in S_merge.free_symbols:
            # Match each symbol with a uniformly random value between 1 and 10; add to dictionary.
            subs_dict[key] = np.random.uniform(1, 10)
        # Using the dictionary, substitute the new values to each symbol
        subs_matrix = S_merge.subs(subs_dict).evalf()
        # Convert SymPy matrix to a NumPy matrix
        subs_matrix_np = np.array(subs_matrix.tolist(), dtype=float)
        # Calculate rank of the new matrix and append to the list
        rank_list.append(np.linalg.matrix_rank(subs_matrix_np))
    return rank_list


# --------------------------------------------------------------
# Create sparse matrix for faster rank solving function for MCM 
# --------------------------------------------------------------
"""
For very large mechanisms like MCM: 
Convert a large symbolic SymPy sparse matrix into a numeric SciPy Compressed Sparse Column (CSC) format matrix.

Args:
    S_merge (SymPy matrix): stoichiometric matrix after merging coproducing reactions

Returns:
    csc_sparse (SciPy sparse matrix): SciPy CSC sparse matrix replica of S_merge 
"""
from scipy.sparse import csc_matrix

def numeric_sparse_matrix_fast_combined(S_merge, seed=None):

    # Optional: set random seed so results are reproducible
    if seed is not None:
        np.random.seed(seed)

    print("Precompute DOK entries")
    # Convert matrix to DOK (Dictionary Of Keys) format to access nonzero entries as (row, col) -> value
    dok = S_merge.todok()
    # Extract positions (i, j) of nonzero entries then convert them into a list
    keys = list(dok.keys())         
    # Extract the list of corresponding symbolic expressions
    values = list(dok.values())  

    # If matrix has no nonzero entries, return an empty sparse matrix
    if not values:
        return csc_matrix(S_merge.shape)

    print("Get all free symbols in order")
    # Collect all symbolic variables appearing anywhere in the matrix and sort them for consistent ordering for function inputs
    free_syms = sorted(list(S_merge.free_symbols), key=lambda s: s.name)

    print("Create fast numeric function using lambdify")
    # Convert symbolic expressions into a fast numerical function
    f = sp.lambdify(free_syms, values, modules='numpy')

    print("Sample random values for symbols")
    # Generate random values for each symbol (uniform between 1 and 10)
    vals = np.random.uniform(1.0, 10.0, size=len(free_syms))

    print("Evaluate all nonzero entries at once")
    # Evaluate all symbolic expressions, creating a numeric array corresponding to each nonzero entry
    data = np.array(f(*vals), dtype=np.float64)

    print("Extract row and column indices")
    # Separate row and column indices from keys
    rows, cols = zip(*keys)

    print("Build CSC sparse matrix")
    # Construct a numeric sparse matrix in CSC format using:
    # - data: evaluated values
    # - (rows, cols): positions of those values
    csc_sparse = csc_matrix((data, (rows, cols)), shape=S_merge.shape)
    return csc_sparse

# ---------------------------------
# Linear algebra function for MCM 
# ---------------------------------
"""
For very large mechanisms like MCM: 
Perform a projection to shrink the matrix to a smaller dimension before computing the rank. 
Repeat for desired number of experiments.

Args:
    A (SciPy sparse matrix): SciPy CSC sparse matrix replica of S_merge; result of calling numeric_sparse_matrix_fast_combined(S_merge, seed=None)
    num_experiments (int): the number of desired experiments
    proj_dim (int): the dimension we want to shrink the matrix to

Returns:
    rank_list (list): a list containing the ranks of the matrix with substituted values, repeated per experiment
"""
def linalg_experiment_fast(A, num_experiments, proj_dim):
    # Get shape of the input matrix
    n_rows, n_cols = A.shape
    # Initialize an empty list
    rank_list = []

    for _ in range(num_experiments):
        # Generate a random projection matrix with size (n_cols, proj_dim)
        R = np.random.randn(n_cols, proj_dim)
        # Project A into the lower-dimensional space while approximately preserving rank structure
        AR = A @ R  

        # Compute the rank of the projected matrix
        rank = np.linalg.matrix_rank(AR, tol=1e-8)
        # Append the rank value to the list
        rank_list.append(rank)

    return rank_list

# ----------------------------------------------------------
# Kinetic Invariant String Constructor (Superfast, Logan81)
# ----------------------------------------------------------
"""
Construct the kinetic invariant vector as a string combination of the involved species, with their respective coefficients.
Used in Superfast and Logan81 mechanisms only.

Args:
    null_vector (SymPy sparse matrix): a sparse matrix (vector) corresponding to the left null space of a merged matrix, obtained within thei ndividual mechanism files
    species_index (dictionary): a species index matching dictionary, created within the individual mechanism files

Returns:
    equation (string): the kinetic invariant formatted as a combination of the involved species
"""
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

