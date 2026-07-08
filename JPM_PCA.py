#%%
import torch
import numpy as np
import pandas as pd

def load_data(folder = './', file = "1000experiments_JPMv1.2_varyingtemperature.csv"):
    
    n_points_per_experiment = 13 # 1 hours * 1 step/5 minutes + 1 for the first step
    n_steps_per_experiment = n_points_per_experiment - 1
    
    # Load reference model output:
    #   C -- concentration values
    
    # Additionally, calculate and output
    #   D -- tendency values (units of concentration)
    print("loading concentration data")
    # read dataframe from csv
    df = pd.read_csv(folder + file,  dtype='float64')

    # Get concentration values C, which are all columns except the first two
    C = df.iloc[:,2:].values

    # Get tendency values D, which are the difference between consecutive concentration values
    D = np.delete(np.diff(C,axis = 0), list(range(n_steps_per_experiment, C.shape[0]-1, n_points_per_experiment)), axis=0)

    # Get C of active species, ignoring H2O, O2, and buildup HNO3, CO, H2
    ignore_species = ['H2O', 'O2', 'HNO3', 'CO', 'H2']
    active_species_columns = [col for col in df.columns if col.split(' ')[0] not in ignore_species and col.endswith('[ppb]')]
    C_active = df[active_species_columns].values

    return df,C,D,C_active

def createIO(C,D,C_active):
    
    n_points_per_experiment = 13 # 1 hours * 1 step/5 minutes + 1 for the first step
    n_steps_per_experiment = n_points_per_experiment - 1
    
    print("creating input and output data")
    X = C_active   
    # delete the last step of each experiment -> Because the last step is the first step of the next experiment! TK
    X = np.delete(C_active, 
                  list(range(n_steps_per_experiment, 
                             C_active.shape[0], n_points_per_experiment)), 
                             axis=0)

    X_all = np.delete(C, 
                  list(range(n_steps_per_experiment, 
                             C_active.shape[0], n_points_per_experiment)), 
                             axis=0)

    Y = D

    #12 consecutive steps are 1 experiment, then the a new series is started!indeces 11,23,35...
    idx_withoutExp_shift = np.array(range(n_steps_per_experiment-1, X.shape[0], n_steps_per_experiment))
    diff = X_all[:-1]-X_all[1:]+D[:-1]
    diff = np.delete(diff, idx_withoutExp_shift[:-1], axis=0)
    print("max diff all:", np.max(np.abs(diff)))
    print("mean diff all:", np.mean(np.abs(diff)))

    # Create a train/test split
    split = 0.90
    trainsplit = int(split*X.shape[0])
    print("train size:", trainsplit)
    testsplit = int(round(1-split,2)*X.shape[0])
    print("test size:", testsplit)
    num_test_exps = int(testsplit/(n_steps_per_experiment))
    print("number of test experiments:", num_test_exps)

    
    X_train_raw = X[0:trainsplit,:]
    Y_train_raw = Y[0:trainsplit,:]
    X_test_raw = X[trainsplit:,:]
    Y_test_raw = Y[trainsplit:,:]
    C_test_raw = C[int(trainsplit*13/12):,:]
    C_test = np.reshape(C_test_raw, [num_test_exps, n_points_per_experiment, C.shape[1]])

    #to torch tensors
    X_train_raw = torch.tensor(X_train_raw)
    Y_train_raw = torch.tensor(Y_train_raw)
    X_test_raw = torch.tensor(X_test_raw)
    Y_test_raw = torch.tensor(Y_test_raw)

    return X_train_raw, Y_train_raw, X_test_raw, Y_test_raw, C_test


#%%
import numpy as np
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.gridspec as gridspec 


# shape: 16 x 13
Svv_JPM3_hardcoded = sp.Matrix([
               [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # O₃
               [ 1, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  1,  0,  0,  0,  1, -1,  0,  0,  0,  0, -1,  1,  1,  0], # NO₂
               [ 0,  0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # HCHO 
               [ 0,  0,  2,  0,  1, -1,  0,  0,  1,  0,  1,  0,  0,  0,  1], # HO₂  
               [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,  0,  0, -1,  1, -1,  2, -1, -1,  0,  0,  0, -1, -1], # OH
               [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0], # CO   
               [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # ALD2
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0], # MGLY
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1,  0,  0], # MCO₃
               [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  0], # PAN
               [ 0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0], # H₂O
               [-1,  1,  0,  0, -1,  0,  0,  0,  0, -1, -2,  0,  0,  0,  1]])# O₂  

a = sp.Symbol('a')

# merge reactions 3 and 4 because of proportional reactants
# shape: 16 x 10, species by reactions
# treat O2 as non existent
S_merge = sp.Matrix([
               [ 1, -1,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1], # O₃
               [ 1, -1,   0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # NO
               [-1,  1,   0,  0,  1, -1,  0,  0,  0,  0, -1,  1,  1,  0], # NO₂
               [ 0,  0,  -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # HCHO 
               [ 0,  0, 2*a,  1, -1,  0,  0,  1,  0,  1,  0,  0,  0,  1], # HO₂  
               [ 0,  0,   0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0], # H₂O₂
               [ 0,  0,   0, -1,  1, -1,  2, -1, -1,  0,  0,  0, -1, -1], # OH
               [ 0,  0,   0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0], # HNO₃
               [ 0,  0,   1,  1,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0], # CO   
               [ 0,  0, 1-a,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # H₂   
               [ 0,  0,   0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # ALD2
               [ 0,  0,   0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0], # MGLY
               [ 0,  0,   0,  0,  0,  0,  0,  0,  1,  1, -1,  1,  0,  0], # MCO₃
               [ 0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  0], # PAN
               [ 0,  0,   0,  1,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0], # H₂O
               [ 0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]])# O₂  

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
LC = M_atom[:, 0]
LN = M_atom[:, 1]
LH = M_atom[:, 2]
LO = M_atom[:, 3]

# check nullspace of merged
L_matrix = sp.Matrix.hstack(*S_merge.T.nullspace())

# append conservation laws to merged matrix
# obtain nullspace
S_merge_conserved = S_merge.row_join(M_atom)
N = sp.Matrix(S_merge_conserved.T.nullspace()[0])

# Combine and simplify all expressions
M_simplified = N.applyfunc(sp.simplify)
# just to check
S_merge_new = S_merge.subs({a:0.598546359448017})
N_new = N.subs({a:0.598546359448017})
S_merge_new.T@N_new

# Multiply through by the least common multiple of denominators to clear fractions (optional)
denominators = [sp.denom(x) for x in M_simplified if not x.is_number]
lcm = sp.lcm(denominators)
M_clean = (M_simplified * lcm).applyfunc(sp.simplify)

simplified_check = (S_merge_conserved.T * M_clean).applyfunc(sp.simplify)
print(simplified_check)

# kinetic invariant
L_alpha_old = sp.Matrix([[a, -a, 0, 0,  a, 0, 0, 0, 0, -2*(1-a),  -a,       a,       0,       0, -a, 0]])
L_alpha =     sp.Matrix([[a, -a, 0, 0,  a, 0, 0, 0, 0, -2*(1-a), a-2, (3*a)-2, (2*a)-2, (2*a)-2, -a, 0]])


# append new kinetic invariant to original conservation laws
M_L_alpha = M_atom.row_join(L_alpha.T)
# check if full rank (linear independence)
print(M_L_alpha.rank())
M_L_alpha_sub = M_L_alpha.subs({a:0.598546359448017})
np_M_L = np.array(M_L_alpha_sub).astype(np.float64)
print(np.linalg.matrix_rank(np_M_L))


#%%

# Check PCA on D data
# SEED = 80
# np.random.seed(SEED)
df, C, D, C_active = load_data(file= '1000experiments_JPMv1.2_varyingtemperature.csv') 
x_train,y_train, x_test, y_test, C_test = createIO(C,D,C_active)
D = D[:, 1:]

# no reduction to check which components vanish
pca2 = PCA(n_components=16, svd_solver="full", random_state=18)   
X_pca_2 = pca2.fit_transform(D)
pca2.singular_values_
y = pca2.singular_values_
y_fixed = np.where(y == 0, 1e-20, y)
explained_variance_ratio = pca2.explained_variance_ratio_

#%%
# Analytical method: Absolute Atom Deviation 

n = len(D)

abs_atm_dev_analytical = D@M_L_alpha_sub
abs_atm_dev_analytical = np.array(abs_atm_dev_analytical, dtype=float)
abs_atm_dev_datadriven = D@pca2.components_[11:, :].T
abs_atm_dev_datadriven = np.array(abs_atm_dev_datadriven, dtype=float)

# Data-driven method: Absolute Atom Deviation

labels = [r'$\Delta [C]$', r'$\Delta [N]$', r'$\Delta [H]$', r'$\Delta [O]$', r'$KI$']
colors = ["green", "#8B4513", "blue", "grey", "#df00e3"]
labels1 = [r'$PC_{11}$', r'$PC_{12}$', r'$PC_{13}$', r'$PC_{14}$', r'$PC_{15}$']
colors1 = ["#34ebde", "#5b00b5", "red", "#cfc800", "#a4ff7a"]


figA = plt.figure(figsize=(8, 11)) 
gs = gridspec.GridSpec(12, 1, figure=figA) 
ax1 = figA.add_subplot(gs[0:6]) 
ax2 = figA.add_subplot(gs[6:9])
ax3 = figA.add_subplot(gs[9:]) 

# Plot PCA Singular Values
c = ["#3f6ea1"]*11
c = c+colors1
ax1.scatter(range(1, 17), y_fixed, c=c, alpha=1.0)
ax1.plot(range(1, 17), y_fixed, color='gray', linewidth=1) 
ax1.axhline(y=1e-5, color='red', linestyle='--', linewidth=1) 
ax1.set_yscale('log') 
ax1.set_xlabel("Principal Components", fontsize=16) 
ax1.set_ylabel(r'$\sigma_i$ : PC Singular Value', fontsize=16) 
ax1.set_xticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
ax1.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

# Plot Atom Deviation of data-driven method
for i in range(abs_atm_dev_analytical.shape[1]): 
    sns.kdeplot(np.log10(abs_atm_dev_datadriven[:, i]), ax=ax2, label=labels1[i], fill=True, alpha=0.5, color=colors1[i], bw_adjust=5) 
    
ax2.set_xlabel("Absolute atom deviation [ppb]", fontsize=16) 
ax2.set_ylabel("Data-Driven Invariants", rotation=90, fontsize=16)
ax2.spines['top'].set_visible(False) 
ax2.spines['right'].set_visible(False) 
ax2.set_yticks([]) 
ax2.set_xticks(np.log10([1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6])) 
ax2.set_xticklabels([r"$10^{-16}$", r"$10^{-14}$", r"$10^{-12}$", r"$10^{-10}$", r"$10^{-8}$", r"$10^{-6}$"]) 
ax2.legend() 
ax2.set_xlim(np.log10(1e-17), np.log10(1e-4))

# Plot Atom Deviation of analytical method
for i in range(abs_atm_dev_analytical.shape[1]): 
    sns.kdeplot(np.log10(abs_atm_dev_analytical[:, i]), ax=ax3, label=labels[i], fill=True, alpha=0.5, color=colors[i], bw_adjust=5) 

ax3.set_xlabel("Absolute atom deviation [ppb]", fontsize=16) 
ax3.set_ylabel("Analytical Invariants", rotation=90, fontsize=16) 
ax3.spines['top'].set_visible(False) 
ax3.spines['right'].set_visible(False) 
ax3.set_yticks([]) 
ax3.set_xticks(np.log10([1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6])) 
ax3.set_xticklabels([r"$10^{-16}$", r"$10^{-14}$", r"$10^{-12}$", r"$10^{-10}$", r"$10^{-8}$", r"$10^{-6}$"])
ax3.legend() 
ax3.set_xlim(np.log10(1e-17), np.log10(1e-4))

ax2.legend(fontsize=11, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)
ax3.legend(fontsize=11, markerscale=0.8, handletextpad=0.4, labelspacing=0.3, frameon=False)

ax1.tick_params(axis='both', labelsize=13)
ax2.tick_params(axis='both', labelsize=13)
ax3.tick_params(axis='both', labelsize=13)

ax1.text(-0.1, 1.05, "(b)", transform=ax1.transAxes,
         fontsize=22, fontweight='bold', va='top', ha='right')

ax2.text(-0.1, 1.05, "(c)", transform=ax2.transAxes,
         fontsize=22, fontweight='bold', va='top', ha='right')

ax3.text(-0.1, 1.05, "(d)", transform=ax3.transAxes,
         fontsize=22, fontweight='bold', va='top', ha='right')

ax1.text(0.7, 0.00005, r'Invariant Threshold ($10^{-5}$)', fontsize=12, color="red")



plt.tight_layout() 
plt.savefig("singular_value_&_atom_deviations.png", dpi=1200)
plt.savefig("singular_value_&_atom_deviations.pdf", dpi=1200)
plt.show()

# %%
