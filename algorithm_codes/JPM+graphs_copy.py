import numpy as np
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from reduce_data import load_data, createIO
import seaborn as sns
import matplotlib.gridspec as gridspec 


# shape: 16 x 13
Svv = sp.Matrix([
               [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O₃
               [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1], # NO₂
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0], # HCHO 
               [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0], # HO₂  
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0], # OH
               [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0], # CO   
               [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1], # PAN
               [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0], # H₂O
               [-1,  1, -2,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0]])# O₂  

# remove reversible reaction redundancies - rxn 2 and 13
# shape: 16 x 11
So = sp.Matrix([
               [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O₃
               [ 1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0], # NO
               [-1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1], # NO₂
               [ 0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0], # HCHO 
               [ 0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0], # HO₂  
               [ 0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0], # H₂O₂
               [ 0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0], # OH
               [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0], # CO   
               [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
               [ 0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1], # PAN
               [ 0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0], # H₂O
               [-1, -2,  0, -1,  0,  0,  0,  0, -1, -2,  0]])# O₂

a = sp.Symbol('a')

# merge reactions 3 and 4 (now w col index 2,3) because of proportional reactants
# shape: 16 x 10, species by reactions
S_merge = sp.Matrix([
                    [ 1,    0,  0,  0,  0,  0,  0,  0,  0,  0], # O₃
                    [ 1,    0,  0, -1,  0,  0,  0,  0,  0,  0], # NO
                    [-1,    0,  0,  1, -1,  0,  0,  0,  0, -1], # NO₂
                    [ 0,   -1, -1,  0,  0,  0,  0,  0,  0,  0], # HCHO 
                    [ 0,  2*a,  1, -1,  0,  0,  1,  0,  1,  0], # HO₂  
                    [ 0,    0,  0,  0,  0, -1, -1,  0,  0,  0], # H₂O₂
                    [ 0,    0, -1,  1, -1,  2, -1, -1,  0,  0], # OH
                    [ 0,    0,  0,  0,  1,  0,  0,  0,  0,  0], # HNO₃
                    [ 0,    1,  1,  0,  0,  0,  0,  0,  1,  0], # CO   
                    [ 0,  1-a,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
                    [ 0,    0,  0,  0,  0,  0,  0, -1,  0,  0], # ALD2
                    [ 0,    0,  0,  0,  0,  0,  0,  0, -1,  0], # MGLY
                    [ 0,    0,  0,  0,  0,  0,  0,  1,  1, -1], # MCO₃
                    [ 0,    0,  0,  0,  0,  0,  0,  0,  0,  1], # PAN
                    [ 0,    0,  1,  0,  0,  0,  1,  1,  0,  0], # H₂O
                    [-1, -2*a, -1,  0,  0,  0,  0, -1, -2,  0]])# O₂  

# merged Rm vector
Rm = [0, 0, 0, -1, 2*a, 0, 0, 0, 1, 1-a, 0, -2, 0, 0, 0, 0]

# treat O2 as non existent
# shape: 16 x 10, species by reactions
S_merge_noO2 = sp.Matrix([
                    [ 1,   0,  0,  0,  0,  0,  0,  0,  0,  0], # O₃
                    [ 1,   0,  0, -1,  0,  0,  0,  0,  0,  0], # NO
                    [-1,   0,  0,  1, -1,  0,  0,  0,  0, -1], # NO₂
                    [ 0,  -1, -1,  0,  0,  0,  0,  0,  0,  0], # HCHO 
                    [ 0, 2*a,  1, -1,  0,  0,  1,  0,  1,  0], # HO₂
                    [ 0,   0,  0,  0,  0, -1, -1,  0,  0,  0], # H₂O₂
                    [ 0,   0, -1,  1, -1,  2, -1, -1,  0,  0], # OH
                    [ 0,   0,  0,  0,  1,  0,  0,  0,  0,  0], # HNO₃
                    [ 0,   1,  1,  0,  0,  0,  0,  0,  1,  0], # CO 
                    [ 0, 1-a,  0,  0,  0,  0,  0,  0,  0,  0], # H₂ 
                    [ 0,   0,  0,  0,  0,  0,  0, -1,  0,  0], # ALD2
                    [ 0,   0,  0,  0,  0,  0,  0,  0, -1,  0], # MGLY
                    [ 0,   0,  0,  0,  0,  0,  0,  1,  1, -1], # MCO₃
                    [ 0,   0,  0,  0,  0,  0,  0,  0,  0,  1], # PAN
                    [ 0,   0,  1,  0,  0,  0,  1,  1,  0,  0], # H₂O
                    [ 0,   0,  0,  0,  0,  0,  0,  0,  0,  0]  # O₂ 
                    ])

# Atom matrix -- SAME AS 4 CONSERVATION LAWS
M_atom = sp.Matrix(
    [
        [0, 0, 0, 3],  # O₃
        [0, 1, 0, 1],  # NO
        [0, 1, 0, 2],  # NO₂
        [1, 0, 2, 1],  # HCHO
        [0, 0, 1, 2],  # HO₂
        [0, 0, 2, 2],  # H₂O₂
        [0, 0, 1, 1],  # OH
        [0, 1, 1, 3],  # HNO₃
        [1, 0, 0, 1],  # CO
        [0, 0, 2, 0],  # H₂
        [2, 0, 4, 1],  # ALD2
        [3, 0, 4, 2],  # MGLY
        [2, 0, 3, 3],  # MCO₃
        [2, 1, 3, 5],  # PAN
        [0, 0, 2, 1],  # H₂O
        [0, 0, 0, 2]   # O₂
    ]
)

# Known conservation laws
LN = M_atom[:, 1]
LH = M_atom[:, 2]
LC = M_atom[:, 0]
LO = M_atom[:, 3]


# check nullspace of merged
L_matrix = sp.Matrix.hstack(*S_merge.T.nullspace())


# append conservation laws to merged matrix
# obtain nullspace
S_merge_conserved = S_merge.row_join(M_atom)
em1 = S_merge_conserved.T.nullspace()[0]
em2 = S_merge_conserved.T.nullspace()[1]
M = sp.Matrix.hstack(em1, em2)


# Combine and simplify all expressions
M_simplified = M.applyfunc(sp.simplify)

# If there are still denominators, bring each element to a common denominator
M_simplified = M_simplified.applyfunc(sp.together)

# Multiply through by the least common multiple of denominators to clear fractions (optional)
denominators = [sp.denom(x) for x in M_simplified if not x.is_number]
lcm = sp.lcm(denominators)
M_clean = (M_simplified * lcm).applyfunc(sp.simplify)

sp.pprint(M_clean)

simplified_check = (S_merge_conserved.T * M_clean).applyfunc(sp.simplify)
print(simplified_check)

# Obin found some invariants
hcho_mgly_co = sp.Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])

# add hcho_mgly_co to S_merge_conserved
# new dimensions: 16x15 --> extra space for the fifth invariant
S_merge_conserve_n5 = S_merge_conserved.row_join(hcho_mgly_co.T)

N = sp.Matrix.hstack(*S_merge_conserve_n5.T.nullspace())
S_merge_new = S_merge.subs({a: 0.5})
N_new = N.subs({a: 0.5})
S_merge_new.T@N_new

# kinetics invariant
L_alpha = sp.Matrix([[a, -a, 0, 0, a, 0, 0, 0, 0, -2*(1-a), -a, a, 0, 0, -a, 0]])

# append new kinetics invariant to original conservation laws
M_atom2 = M_atom.row_join(hcho_mgly_co.T)
M_L_alpha = M_atom2.row_join(L_alpha.T)
# check if full rank (linear independence)
print(M_L_alpha.rank())
M_L_alpha_sub = M_L_alpha.subs({a:0.598546359448017})
np_M_L = np.array(M_L_alpha_sub).astype(np.float64)
print(np.linalg.matrix_rank(np_M_L))

# Check PCA on D data
SEED = 42
np.random.seed(SEED)

df, C, D, C_active = load_data(file= 'experiments_11e5_1hour_5mins_falsecombinatoricratelaws.csv') #replace with your path
x_train,y_train, x_test, y_test, C_test = createIO(C,D,C_active)

# no reduction to check which components vanish
pca = PCA(n_components=16, svd_solver="full", random_state=18)   
X_pca = pca.fit_transform(D)
pca.singular_values_
y = pca.singular_values_
y_fixed = np.where(y == 0, 1e-20, y)
explained_variance_ratio = pca.explained_variance_ratio_

# Analytical method: Absolute Atom Deviation 

abs_atm_dev_analytical = D[:50000, :]@M_L_alpha_sub
abs_atm_dev_analytical = np.array(abs_atm_dev_analytical, dtype=float)
abs_atm_dev_datadriven = D[:50000, :]@pca.components_[10:, :].T
abs_atm_dev_datadriven = np.array(abs_atm_dev_datadriven, dtype=float)

# Data-driven method: Absolute Atom Deviation

labels = ["C Conservation", "N Conservation", "H Conservation", "O Conservation", "C1 Carbon Subpool Invariant", "Coproduction Emanant"]
colors = ["green", "#8B4513", "blue", "grey", "orange", "purple"]
labels1 = ["PC11", "PC12", "PC13", "PC14", "PC15", "PC16"]
colors1 = ["#8ccdff", "#cfc800", "red", "#7251a6", "#df00e3", "#00a881"]

figA = plt.figure(figsize=(7, 9)) 
gs = gridspec.GridSpec(8, 1, figure=figA) 
ax1 = figA.add_subplot(gs[0:4]) 
ax2 = figA.add_subplot(gs[4:6])
ax3 = figA.add_subplot(gs[6:]) 

# Plot PCA Singular Values
c = ["#3f6ea1"]*10
c = c+colors1
ax1.scatter(range(1, 17), y_fixed, c=c, alpha=1.0)
ax1.plot(range(1, 17), y_fixed, color='gray', linewidth=1) 
ax1.axhline(y=1e-6, color='red', linestyle='--', linewidth=1) 
ax1.set_yscale('log') 
ax1.set_xlabel("Principal Components") 
ax1.set_ylabel(r'$\sigma_i$ : PC Singular Value') 
ax1.set_xticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
ax1.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])

# Plot Atom Deviation of data-driven method
for i in range(abs_atm_dev_analytical.shape[1]): 
    sns.kdeplot(np.log10(abs_atm_dev_datadriven[:, i]), ax=ax2, label=labels1[i], fill=True, alpha=0.5, color=colors1[i], bw_adjust=5) 
    
ax2.set_xlabel("Absolute atom deviation [ppb]") 
ax2.set_ylabel("Data-Driven Invariants", rotation=90)
ax2.spines['top'].set_visible(False) 
ax2.spines['right'].set_visible(False) 
ax2.set_yticks([]) 
ax2.set_xticks(np.log10([1e-20, 1e-15, 1e-10, 1e-5, 1])) 
ax2.set_xticklabels([r"$10^{-20}$", r"$10^{-15}$", r"$10^{-10}$", r"$10^{-5}$", r"$10^{0}$"]) 
ax2.legend() 

# Plot Atom Deviation of analytical method
for i in range(abs_atm_dev_analytical.shape[1]): 
    sns.kdeplot(np.log10(abs_atm_dev_analytical[:, i]), ax=ax3, label=labels[i], fill=True, alpha=0.5, color=colors[i], bw_adjust=5) 

ax3.set_xlabel("Absolute atom deviation [ppb]") 
ax3.set_ylabel("Analytical Invariants", rotation=90) 
ax3.spines['top'].set_visible(False) 
ax3.spines['right'].set_visible(False) 
ax3.set_yticks([]) 
ax3.set_xticks(np.log10([1e-20, 1e-15, 1e-10, 1e-5, 1])) 
ax3.set_xticklabels([r"$10^{-20}$", r"$10^{-15}$", r"$10^{-10}$", r"$10^{-5}$", r"$10^{0}$"])
ax3.legend() 

ax2.legend(fontsize=8, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)
ax3.legend(fontsize=8, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)

ax1.text(-0.1, 1.05, "(b)", transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')

ax2.text(-0.1, 1.05, "(c)", transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')

ax3.text(-0.1, 1.05, "(d)", transform=ax3.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')


plt.tight_layout() 
plt.show()
plt.savefig("singular_value_&_atom_deviations.png")
