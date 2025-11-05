import numpy as np
import sympy as sp
from scipy.linalg import null_space


# Species order:
species = [
    "O3", "NO", "NO2", "HCHO", "HO2", "H2O2", "OH", "HNO3",
    "CO", "H2", "H2O", "O2", "ALD2", "MGLY", "MCO3", "PAN"
]

# Reaction order:
reactions = [
    "rk1", "rk2", "rk3", "rk4", "rk5", "rk6", "rk7",
    "rk8", "rk9", "rk10", "rk11", "rk12", "rk13"
]

# Stoichiometric matrix S (rows = species, columns = reactions)
S = np.array([
    [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # O3
    [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],  # NO
    [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1],  # NO2
    [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0],  # HCHO
    [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0],  # HO2
    [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0],  # H2O2
    [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0],  # OH
    [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],  # HNO3
    [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0],  # CO
    [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # H2
    [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0],  # H2O
    [-1,  1, -2,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0],  # O2
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],  # ALD2
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],  # MGLY
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1],  # MCO3
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1],  # PAN
], dtype=float)


# Atom order: [C, N, H, O]
atoms = np.array([
    [0, 0, 0, 3],  # O3
    [0, 1, 0, 1],  # NO
    [0, 1, 0, 2],  # NO2
    [1, 0, 2, 1],  # HCHO
    [0, 0, 1, 2],  # HO2
    [0, 0, 2, 2],  # H2O2
    [0, 0, 1, 1],  # OH
    [0, 1, 1, 3],  # HNO3
    [1, 0, 0, 1],  # CO
    [0, 0, 2, 0],  # H2
    [0, 0, 2, 1],  # H2O
    [0, 0, 0, 2],  # O2
    [2, 0, 4, 1],  # ALD2
    [3, 0, 4, 2],  # MGLY
    [2, 0, 3, 3],  # MCO3
    [2, 1, 3, 5],  # PAN
], dtype=float)

hcho_mgly_co = np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0])
ald2_mgly_mco3_pan = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])

# Quick check of conservation
conservation = atoms.T @ S  # Shape (4, number of reactions)
# check rank of [atoms, hcho_mgly_co]
AH = np.hstack((atoms, np.array(hcho_mgly_co).reshape(-1, 1)))
sp.Matrix(AH).rank() 

# potential conserved quantities
ns_pot = np.array([-21,  -610,  -617,  -5280,  -460,  -906,  -453,  -1070,
      -4388, -892,  -899,  -14,   2482, -1906,  2914,   2297])
ns_pot.T@S  # should be close to zero

nullspace = sp.Matrix.nullspace(sp.Matrix(S.T))
nullspace = sp.Matrix.hstack(*nullspace)


# what is in the nullspace of [S, atoms]?
NS = sp.Matrix.nullspace(sp.Matrix(np.vstack((S.T, atoms.T))))
NS = sp.Matrix.hstack(*NS)
# get the least common multiple of denominators, within tolerance of .0001
# lcm_denoms = sp.ilcm(*[v.q for v in NS])
# # scale to integer vectors
# NS_int = NS.applyfunc(lambda x: int(x * lcm_denoms))

invariants = np.hstack([atoms, ns_pot.reshape(-1, 1)])  # shape (4, 17)

# solve nullspace*X = invariants (or atoms)
# do this with a right inverse
nullspace_np = np.array(nullspace).astype(np.float64)
nullspace_pinv = np.linalg.pinv(nullspace_np)
X = nullspace_pinv @ atoms.astype(np.float64)


nullspace_X = sp.Matrix.nullspace(sp.Matrix(X.T))[0]
nullspace_X = nullspace_X/3.6e16*4 

# get the 5th emergent conservation law
emergent = nullspace@nullspace_X

# round emergent to integers (nearest 1)
emergent_rounded = np.array(emergent).astype(np.float64)
emergent_rounded = np.round(emergent_rounded).astype(int)

# what's the rank of [atoms, emergent]?
AE = np.hstack((atoms, emergent_rounded.reshape(-1, 1)))
sp.Matrix(AE).rank()  # should be 5