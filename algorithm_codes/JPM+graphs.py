import numpy as np
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
from reduce_data import load_data, createIO
import seaborn as sns
import matplotlib.gridspec as gridspec 


# shape: 16 x 13
Svv = sp.Matrix([
               [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [-1, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 1],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, -1, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 1, -1, 2, -1, -1, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
               [-1, 1, -2, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0]])

# remove reversible reaction redundancies - rxn 2 and 13
# shape: 16 x 11
So = sp.Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1],
               [0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
               [0, 2, 0, 1, -1, 0, 0, 1, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0],
               [0, 0, 0, -1, 1, -1, 2, -1, -1, 0, 0],
               [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
               [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0],
               [-1, -2, 0, -1, 0, 0, 0, 0, -1, -2, 0]])

a = sp.Symbol('a')

# merge reactions 3 and 4 (now w col index 2,3) because of proportional reactants
# shape: 16 x 10, species by reactions
S_merge = sp.Matrix([
                    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],       # O₃
                    [1, 0, 0, -1, 0, 0, 0, 0, 0, 0],      # NO
                    [-1, 0, 0, 1, -1, 0, 0, 0, 0, -1],    # NO₂
                    [0, -1, -1, 0, 0, 0, 0, 0, 0, 0],     # HCHO --> 3
                    [0, 2*a, 1, -1, 0, 0, 1, 0, 1, 0],    # HO₂  --> 4
                    [0, 0, 0, 0, 0, -1, -1, 0, 0, 0],     # H₂O₂
                    [0, 0, -1, 1, -1, 2, -1, -1, 0, 0],   # OH
                    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],       # HNO₃
                    [0, 1, 1, 0, 0, 0, 0, 0, 1, 0],       # CO   --> 8
                    [0, 1-a, 0, 0, 0, 0, 0, 0, 0, 0],     # H₂   --> 9
                    [0, 0, 0, 0, 0, 0, 0, -1, 0, 0],      # ALD2
                    [0, 0, 0, 0, 0, 0, 0, 0, -1, 0],      # MGLY
                    [0, 0, 0, 0, 0, 0, 0, 1, 1, -1],      # MCO₃
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],       # PAN
                    [0, 0, 1, 0, 0, 0, 1, 1, 0, 0],       # H₂O
                    [-1, -2*a, -1, 0, 0, 0, 0, -1, -2, 0],# O₂   --> 15
                    ])

# Left nullspace of S = conserved quantities
# Right nullspace of S.T = conserved quantities
# rank = s - l = r - c
# for Svv:
# 11 = 16 - 5 = 13 - 2
# for So:
# 11 = 16 - 5 = 11 - 0
# for S_merge:
# 10 = 16 - 6 = 10 - 0
# emanants = 6 - 5 = 1

# condensed merged vector Rm: [-1, 2*a, 1, 1-a, -2]
# merged Rm vector
Rm = [0, 0, 0, -1, 2*a, 0, 0, 0, 1, 1-a, 0, -2, 0, 0, 0, 0]


# treat O2 as non existent
# shape: 16 x 10, species by reactions
S_merge_noO2 = sp.Matrix([
                    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],       # O₃
                    [1, 0, 0, -1, 0, 0, 0, 0, 0, 0],      # NO
                    [-1, 0, 0, 1, -1, 0, 0, 0, 0, -1],    # NO₂
                    [0, -1, -1, 0, 0, 0, 0, 0, 0, 0],     # HCHO --> 3
                    [0, 2*a, 1, -1, 0, 0, 1, 0, 1, 0],    # HO₂  --> 4
                    [0, 0, 0, 0, 0, -1, -1, 0, 0, 0],     # H₂O₂
                    [0, 0, -1, 1, -1, 2, -1, -1, 0, 0],   # OH
                    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],       # HNO₃
                    [0, 1, 1, 0, 0, 0, 0, 0, 1, 0],       # CO   --> 8
                    [0, 1-a, 0, 0, 0, 0, 0, 0, 0, 0],     # H₂   --> 9
                    [0, 0, 0, 0, 0, 0, 0, -1, 0, 0],      # ALD2
                    [0, 0, 0, 0, 0, 0, 0, 0, -1, 0],      # MGLY
                    [0, 0, 0, 0, 0, 0, 0, 1, 1, -1],      # MCO₃
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],        # PAN
                    [0, 0, 1, 0, 0, 0, 1, 1, 0, 0],       # H₂O
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]        # O₂   --> 15
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

# check if current conservation laws are a combination of the nullspace
# x = sp.linsolve((L_matrix, K_matrix))
# print(x)

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

# sp.pprint(M_clean)

simplified_check = (S_merge_conserved.T * M_clean).applyfunc(sp.simplify)
# print(simplified_check)

# Obin found some invariants
hcho_mgly_co = sp.Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])
# ald2_mgly_mco3_pan = sp.Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0]])
# carbon_2ndhalf = LC - hcho_mgly_co

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
# print(M_L_alpha.rank())
M_L_alpha_sub = M_L_alpha.subs({a:0.598546359448017})
# np_M_L = np.array(M_L_alpha_sub).astype(np.float64)
# print(np.linalg.matrix_rank(np_M_L))

# check if kinetics invariant is in the nullspace of S_merge
# val_alpha = 0.4
# test_S_merge = S_merge.subs({a:val_alpha})
# L_alpha_new = L_alpha.subs({a:1-val_alpha})
# print(test_S_merge.T@L_alpha_new.T)

# CHANGE "a" to "1-a" in S_merge


# Check PCA on D data
SEED = 42
np.random.seed(SEED)

df, C, D, C_active = load_data(file= 'experiments_11e5_1hour_5mins_falsecombinatoricratelaws.csv') #replace with your path
x_train,y_train, x_test, y_test, C_test = createIO(C,D,C_active)

# check conservation in delC data
# cons_check = D[:1000, :]@M_L_alpha_sub

# D_pandas = pd.DataFrame(D)
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(D_pandas)
# cov_matrix = np.cov(X_scaled.T)  
# print("Covariance matrix:\n", cov_matrix)

# eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
# print("Eigenvalues:\n", eig_vals)
# print("Eigenvectors:\n", eig_vecs)

# idx = np.argsort(eig_vals)[::-1]  
# eig_vals = eig_vals[idx]
# eig_vecs = eig_vecs[:, idx]
# print("Sorted Eigenvalues:\n", eig_vals)

# explained_variance = eig_vals / np.sum(eig_vals)
# cumulative_variance = np.cumsum(explained_variance)
# print("Explained variance per PC:\n", explained_variance)
# print("Cumulative variance:\n", cumulative_variance)

# pca = PCA(n_components=10)   
# X_pca = pca.fit_transform(D)
# pca_df = pd.DataFrame(data=X_pca, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"])

# loadings = pd.DataFrame(
#     pca.components_.T,
#     columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"],
#     index=D_pandas.columns)
# print(loadings)
# print("Component Singular Values:", pca.singular_values_)


# plt.bar(range(1, len(pca.singular_values_)+1), pca.singular_values_)
# plt.xlabel("Principal Component")
# plt.ylabel("Singular Value")
# plt.title("PCA Singular Values")
# plt.show()

# no reduction to check which components vanish
pca2 = PCA(n_components=16, svd_solver="full", random_state=18)   
X_pca_2 = pca2.fit_transform(D)
pca2.singular_values_
y = pca2.singular_values_
y_fixed = np.where(y == 0, 1e-20, y)
explained_variance_ratio = pca2.explained_variance_ratio_
# print(sum(explained_variance_ratio[11:]))

# plt.plot(range(1, 17), y_fixed, marker='o')
# plt.axhline(y=1e-2, color='red', linestyle='--', linewidth=1)
# plt.yscale('log')
# plt.xlabel("Component Index")
# plt.ylabel(r'$\sigma_i$ : Component PCA Singular Value')
# plt.show()
# plt.savefig("PCA_Singular_Values.png")

# def plot_biplot(ax, pca_scores, loadings, labels, pc1_idx, pc2_idx):
#     ax.scatter(pca_scores[:, pc1_idx], pca_scores[:, pc2_idx], alpha=0.5)

#     scaling_factor = X_pca[:, [pc1_idx, pc2_idx]].max() * 0.7
#     scaled_loadings = loadings * scaling_factor

#     for i, label in enumerate(labels):
#         ax.arrow(0, 0, scaled_loadings[i, pc1_idx], scaled_loadings[i, pc2_idx],
#                     color='r', alpha=0.8, head_width=0.05, head_length=0.05)
#         ax.text(scaled_loadings[i, pc1_idx] * 1.1, scaled_loadings[i, pc2_idx] * 1.1,
#                 label, color='r', fontsize=9)

#     ax.set_xlabel(f'Principal Component {pc1_idx + 1}')
#     ax.set_ylabel(f'Principal Component {pc2_idx + 1}')
#     ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
#     ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
#     ax.set_title(f'Biplot: PC{pc1_idx+1} vs PC{pc2_idx+1}')

# fig, axes = plt.subplots(1, 2, figsize=(12, 6))
# labels = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]

# plot_biplot(axes[0], X_pca, loadings.values, labels, 0, 1)
# plot_biplot(axes[1], X_pca, loadings.values, labels, 1, 2)


# Analytical method: Absolute Atom Deviation 

abs_atm_dev_analytical = D[:50000, :]@M_L_alpha_sub
abs_atm_dev_analytical = np.array(abs_atm_dev_analytical, dtype=float)
abs_atm_dev_datadriven = D[:50000, :]@pca2.components_[10:, :].T
abs_atm_dev_datadriven = np.array(abs_atm_dev_datadriven, dtype=float)

# fig, ax = plt.subplots(figsize=(16, 4))

# for i in range(abs_atm_dev_analytical.shape[1]):
#     sns.kdeplot(np.log10(abs_atm_dev_analytical[:, i]), ax=ax, label=labels[i],
#                 fill=True, alpha=0.5, color=colors[i], bw_adjust=5)

# ax.set_xlabel("Absolute atom deviation [ppb]", fontsize=15)
# ax.set_ylabel("Analytical Invariants (Linear Algebra)", rotation=90, fontsize=17, labelpad=12)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.set_yticks([])
# ax.set_xticks(np.log10([1e-20, 1e-15, 1e-10, 1e-5, 1]))
# ax.set_xticklabels([r"$10^{-20}$", r"$10^{-15}$", r"$10^{-10}$", r"$10^{-5}$", r"$10^{0}$"],
#                     fontsize=13)

# plt.legend()
# plt.tight_layout()
# plt.savefig("atom_deviation_analytical.png", dpi=300)
# plt.show()

#%%

# Data-driven method: Absolute Atom Deviation

labels = ["C Conservation", "N Conservation", "H Conservation", "O Conservation", "C1 Carbon Subpool Invariant", "Coproduction Emanant"]
colors = ["green", "#8B4513", "blue", "grey", "orange", "purple"]
labels1 = ["PC11", "PC12", "PC13", "PC14", "PC15", "PC16"]
colors1 = ["#8ccdff", "#cfc800", "red", "#7251a6", "#df00e3", "#00a881"]

# fig1, ax1 = plt.subplots(figsize=(16,4))

# for i in range(abs_atm_dev_analytical.shape[1]):
#     sns.kdeplot(np.log10(abs_atm_dev_datadriven[:, i]), ax=ax1, label=labels1[i],
#                 fill=True, alpha=0.5, color=colors1[i], bw_adjust=5)
    
# ax1.set_xlabel("Absolute atom deviation [ppb]", fontsize=15)
# ax1.set_ylabel("Data-Driven Invariants (PCA)", rotation=90, fontsize=17, labelpad=12)
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.set_yticks([])
# ax1.set_xticks(np.log10([1e-20, 1e-15, 1e-10, 1e-5, 1]))
# ax1.set_xticklabels([r"$10^{-20}$", r"$10^{-15}$", r"$10^{-10}$", r"$10^{-5}$", r"$10^{0}$"],
#                     fontsize=13)

# plt.legend()
# plt.tight_layout()
# plt.savefig("atom_deviation_datadriven.png", dpi=300)
# plt.show()


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


# %%

import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

# I removed the -2 for O2, rxn3 only in Svv and Sr 
Svv = sp.Matrix([
               [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [-1, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 1],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, -1, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 1, -1, 2, -1, -1, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
               [-1, 1, 0, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0]])
Svv_sparse = sp.SparseMatrix(Svv)

Sp = sp.Matrix([
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
Sp_sparse = sp.SparseMatrix(Sp)

Sr = sp.Matrix([
               [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 0, -1, 0, -1, -1, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0]])
Sr_sparse = sp.SparseMatrix(Sr)

stoichiometric_invariants = len(Svv_sparse.T.nullspace())
coproduction_cols = create_coproduction(Sr_sparse)
symbol_dict = create_symbols(coproduction_cols)
S_merge = merge_coprod(Sr_sparse, Sp_sparse, symbol_dict, coproduction_cols, Svv_sparse)

S_merge[15, 2] = -2*symbol_dict[3]

del_l, del_r, del_c = s_linalg(Svv_sparse, S_merge, stoichiometric_invariants)

# %%

import pandas as pd
import numpy as np
import sympy as sp
from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

# add cycle -- reaction 14, coproduces with reaction 7, cycle with reaction 6

# shape 16x14
Svv = sp.Matrix([
               [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  1],
               [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1, -1],
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0,  1],
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0, -1], 
               [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0],
               [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1,  0],
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0],
               [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0,  0],
               [-1,  1,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0,  0]])
Svv_sparse = sp.SparseMatrix(Svv)

# shape 16x14
Sp = sp.Matrix([
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
               [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0], 
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
Sp_sparse = sp.SparseMatrix(Sp)

# shape 16x14
Sr = sp.Matrix([
               [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, -1],
               [0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 0, -1, 0, -1, -1, 0, 0, 0, -1], 
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [-1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -2, 0, 0, 0]])
Sr_sparse = sp.SparseMatrix(Sr)

# 5
stoichiometric_invariants = len(Svv_sparse.T.nullspace())
coproduction_cols = create_coproduction(Sr_sparse)
symbol_dict = create_symbols(coproduction_cols)
S_merge = merge_coprod(Sr_sparse, Sp_sparse, symbol_dict, coproduction_cols, Svv_sparse)

S_merge[15, 2] = -2*symbol_dict[3]

# expect del_r = -2, del_l = 1, del_c = 1
del_l, del_r, del_c = s_linalg(Svv_sparse, S_merge, stoichiometric_invariants)
# %%
