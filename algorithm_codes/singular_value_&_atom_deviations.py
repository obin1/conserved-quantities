import numpy as np
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from algorithm_codes.reduce_data import load_data, createIO
import seaborn as sns
import matplotlib.gridspec as gridspec 


# -------------------------------------------------------
# Defining Stoichiometric Matrices + Merging Procedure
# -------------------------------------------------------

# Define stoichiometric matrix
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

# Remove reversible reaction redundancies - reactions 2 and 13
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

# Define SymPy symbolic variable "a" for the partial kinetic reaction rates
a = sp.Symbol('a')

# merge original reactions 3 and 4 (now with column indices 2 and 3) because of proportional reactants
# shape: 16 x 10
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

# Treat O2 as non existent
# shape: 16 x 10
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


# -------------------------------------------------------
# Conservation Laws (Invariants) Derived Analytically
# -------------------------------------------------------

# Define the atom matrix -- equivalent to 4 stoichiometric invariants in the system
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

# Stoichiometric invariants 
LN = M_atom[:, 1] # Nitrogen conservation
LH = M_atom[:, 2] # Hydrogen conservation
LC = M_atom[:, 0] # Carbon conservation
LO = M_atom[:, 3] # Oxygen conservation

# Discovered Carbon Subpool Invariant (5th stoichiometric invariant) in vector form
hcho_mgly_co = sp.Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])

# Emergent kinetic invariant in vector form
L_alpha = sp.Matrix([[a, -a, 0, 0, a, 0, 0, 0, 0, -2*(1-a), -a, a, 0, 0, -a, 0]])

# Add the carbon subpool invariant and the emergent kinetic invariant to the stoichiometric invariant matrix
M_atom2 = M_atom.row_join(hcho_mgly_co.T)
M_L_alpha = M_atom2.row_join(L_alpha.T)

# Using the reaction rates for R3 and R4, calculate the partial kinetic reaction rate to be 0.598546359448017
# Substitute that value into "a" in the SymPy symbolic matrix
M_L_alpha_sub = M_L_alpha.subs({a:0.598546359448017})


# -------------------------------------------------------
# Principal Component Analysis on Species Tendencies
# -------------------------------------------------------

# Initialize random seed
SEED = 42
np.random.seed(SEED)

# Load simulated time series data of species concentration 
# D: change in species concentration (tendencies) between time points of each experiment
# Structure: 13 consecutive time-stamps (each lasting 5 minutes) are 1 experiment; 
# D contains the 12 tendencies between the 13 time-stamps per experiment
df, C, D, C_active = load_data(file= 'experiments_11e5_1hour_5mins_falsecombinatoricratelaws.csv') 

# Perform PCA on D
pca = PCA(n_components=16, svd_solver="full", random_state=18)   
X_pca = pca.fit_transform(D)
pca.singular_values_
y = pca.singular_values_
y_fixed = np.where(y == 0, 1e-20, y)
explained_variance_ratio = pca.explained_variance_ratio_

# -------------------------------------------------------
# Analytical method: Absolute Atom Deviation 
# -------------------------------------------------------

# Take 50,000 experiments in D, and perform matrix multiplication with the matrix of stoichiometric and kinetic invariants
abs_atm_dev_analytical = D[:50000, :]@M_L_alpha_sub
# Convert the product to a NumPy matrix
abs_atm_dev_analytical = np.array(abs_atm_dev_analytical, dtype=float)

# -------------------------------------------------------
# Data-driven method (PCA): Absolute Atom Deviation
# -------------------------------------------------------

# Take 50,000 experiments in D, and perform matrix multiplication with the vanishing component vectors of the PCA
abs_atm_dev_datadriven = D[:50000, :]@pca.components_[10:, :].T
# Convert the product to a NumPy matrix
abs_atm_dev_datadriven = np.array(abs_atm_dev_datadriven, dtype=float)


# -------------------------------------------------------
# Subplots b, c, d for Figure 2
# -------------------------------------------------------

# Define specific labels and color choices for each invariant
labels = ["C Conservation", "N Conservation", "H Conservation", "O Conservation", "C1 Carbon Subpool Invariant", "Coproduction Emanant"]
colors = ["green", "#8B4513", "blue", "grey", "orange", "purple"]
labels1 = ["PC11", "PC12", "PC13", "PC14", "PC15", "PC16"]
colors1 = ["#8ccdff", "#cfc800", "red", "#7251a6", "#df00e3", "#00a881"]

# Define an empty MatPlotLib figure
figA = plt.figure(figsize=(7, 9)) 
# Create a GridSpec to divide the plot into 3 parts with relative size specifications
gs = gridspec.GridSpec(8, 1, figure=figA) 
ax1 = figA.add_subplot(gs[0:4]) 
ax2 = figA.add_subplot(gs[4:6])
ax3 = figA.add_subplot(gs[6:]) 

# Subplot 1:
# Plot the PCA Singular Values to highlight the existence and the number of vanishing components
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

# Subplot 2:
# Plot the distribution of atomic deviation of the data-driven method
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

# Subplot 3:
# Plot the distribution of atomic deviation of the analytical method
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

# Formatting specification for the legend
ax2.legend(fontsize=8, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)
ax3.legend(fontsize=8, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)

# b), c), and d) labels for the subplots, as seen in Figure 2 of the paper
ax1.text(-0.1, 1.05, "(b)", transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')
ax2.text(-0.1, 1.05, "(c)", transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')
ax3.text(-0.1, 1.05, "(d)", transform=ax3.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='right')

# Display and save image as a png
plt.tight_layout() 
plt.show()
plt.savefig("singular_value_&_atom_deviations.png")
