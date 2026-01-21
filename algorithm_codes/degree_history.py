import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sp
from collections import Counter

from algorithm import create_sparse

edge_list_superfast = pd.read_csv("../mechanisms/superfast/superfast_EdgeList.csv", comment="!")
edge_list_pollu = pd.read_csv("../mechanisms/pollu/pollu_EdgeList.csv", comment="!")
edge_list_saprc99 = pd.read_csv("../mechanisms/saprc99/saprc99_EdgeList.csv", comment="!")
edge_list_small_strato = pd.read_csv("../mechanisms/small_strato/small_strato_EdgeList.csv", comment="!")
edge_list_amore = pd.read_csv("../mechanisms/amore2_isop_135species/amore2_isop_135species_EdgeList.csv", comment="!")
edge_list_cb05 = pd.read_csv("../mechanisms/cb05/CB05TUCl_EPA_EdgeList.csv", comment="!")
edge_list_cbmz = pd.read_csv("../mechanisms/cbmz/CBMZ_EdgeList.csv", comment="!")
edge_list_cri = pd.read_csv("../mechanisms/cri-v2.2/cri_EdgeList.csv", comment="!")
edge_list_chemuci = pd.read_csv("../mechanisms/e3sm-chem/chemuci_EdgeList.csv", comment="!")
edge_list_form1985 = pd.read_csv("../mechanisms/form1985/form1985_EdgeList.csv", comment="!")
edge_list_gc12_0_0 = pd.read_csv("../mechanisms/geoschem-12.0.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc12_7_0 = pd.read_csv("../mechanisms/geoschem-12.7.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc12_8_0 = pd.read_csv("../mechanisms/geoschem-12.8.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc12_9_0 = pd.read_csv("../mechanisms/geoschem-12.9.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc13_3_0 = pd.read_csv("../mechanisms/geoschem-13.3.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc13_4_0 = pd.read_csv("../mechanisms/geoschem-13.4.0/gckpp_EdgeList.csv", comment="!")
edge_list_gc = pd.read_csv("../mechanisms/geos-chem-v14/gckpp_EdgeList.csv", comment="!")
edge_list_gc14_6_3 = pd.read_csv("../mechanisms/geoschem-14.6.3/gckpp_EdgeList.csv", comment="!")
edge_list_gchg= pd.read_csv("../mechanisms/geoschem-hg/Hg_EdgeList.csv", comment="!")
edge_list_iso_full = pd.read_csv("../mechanisms/isoprene_full_v5_nofixedspc/isoprene_full_v5_nofixedspc_EdgeList.csv", comment="!")
edge_list_iso_red = pd.read_csv("../mechanisms/isoprene_reduced_plus_v5/isoprene_reduced_plus_v5_EdgeList.csv", comment="!")
edge_list_mecca = pd.read_csv("../mechanisms/mecca_v4.6.0/gas_EdgeList.csv", comment="!")
edge_list_mcm = pd.read_csv("../mechanisms/mcm_v3.3.1/mcm_EdgeList.csv", comment="!")
edge_list_mozart_t1 = pd.read_csv("../mechanisms/mozart-t1/MOZART_T1_EdgeList.csv", comment="!")
edge_list_mozart = pd.read_csv("../mechanisms/mozart/MOZART_4_EdgeList.csv", comment="!")
edge_list_pact1d = pd.read_csv("../mechanisms/pact1d/mech_EdgeList.csv", comment="!")
edge_list_racm = pd.read_csv("../mechanisms/racm/RACM_EdgeList.csv", comment="!")
edge_list_radm2 = pd.read_csv("../mechanisms/radm2/RADM2_EdgeList.csv", comment="!")

edgelists = [edge_list_superfast,edge_list_pollu,edge_list_saprc99,edge_list_small_strato,edge_list_amore,edge_list_cb05,edge_list_cbmz,edge_list_cri,edge_list_chemuci,edge_list_form1985,edge_list_gc12_0_0,edge_list_gc12_7_0,edge_list_gc12_8_0,edge_list_gc12_9_0,edge_list_gc13_3_0,edge_list_gc13_4_0,edge_list_gc,edge_list_gc14_6_3,edge_list_gchg,edge_list_iso_full,edge_list_iso_red,edge_list_mecca,edge_list_mcm,edge_list_mozart_t1,edge_list_mozart,edge_list_pact1d,edge_list_racm,edge_list_radm2]

mechanisms = ["superfast", "pollu", "saprc99", "small_strato", "amore", "cb05", "cbmz", "cri", "chemuci", "form1985", "gc12.0.0", "gc12.7.0", "gc12.8.0", "gc12.9.0", "gc13.3.0", "gc13.4.0", "gc", "gc14.6.3", "gchg", "iso_full", "iso_red", "mecca", "mcm", "mozart_t1", "mozart", "pact1d", "racm", "radm2"]

#%%

Svv_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O2
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
Svv_sparse_d = sp.SparseMatrix(Svv_d)

Sr_d = sp.Matrix([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O2
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
Sr_sparse_d = sp.SparseMatrix(Sr_d)

Sp_d = sp.Matrix([
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # PhCH2CHO
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # O2
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
Sp_sparse_d = sp.SparseMatrix(Sp_d)

#%%

Svv_JPM = sp.Matrix([
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
Svv_sparse_JPM = sp.SparseMatrix(Svv_JPM)

Sp_JPM = sp.Matrix([
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
Sp_sparse_JPM = sp.SparseMatrix(Sp_JPM)

Sr_JPM = sp.Matrix([
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
Sr_sparse_JPM = sp.SparseMatrix(Sr_JPM)

#%%

edge_dict = {}

for name, e in zip(mechanisms, edgelists):
    Sr, Sp, Svv = create_sparse(e)
    edge_dict[name] = [Sr, Sp, Svv]

edge_dict["daniel"] = [Sr_sparse_d, Sp_sparse_d, Svv_sparse_d]
edge_dict["jpm"] = [Sr_sparse_JPM, Sp_sparse_JPM, Svv_sparse_JPM]

mechanisms.append("daniel")
mechanisms.append("jpm")

count_dict = {}
for name, mech in zip(mechanisms, edge_dict.keys()):
    Sr, Sp, Svv = edge_dict[mech]
    all_list = []
    for i in range(Sr.shape[0]):
        spc_list = []
        for j in range(Sr.shape[1]):
            if Sr[i,j]!=0:
                spc_list.append(j)
            if (Sp[i,j]!=0) & (j not in spc_list):
                spc_list.append(j)
        all_list.append(len(spc_list))
    count_dict[name] = all_list

#%%
plt.figure(figsize=(8, 4))
for mech, counts in count_dict.items():
    freq = Counter(counts)

    x = np.array(sorted(freq.keys()))
    y = np.array([freq[k] for k in x])

    plt.plot(x, y, marker='o', label=mech)

plt.xlabel("Degree / Count")
plt.ylabel("Frequency")
# plt.legend()
plt.tight_layout()
plt.show()


#%%
n_mechs = len(count_dict)
n_cols = 3                  
n_rows = math.ceil(n_mechs / n_cols)

fig, axes = plt.subplots(
    n_rows, n_cols,
    figsize=(4 * n_cols, 3 * n_rows),
)

axes = axes.flatten()
for ax, (mech, counts) in zip(axes, count_dict.items()):
    ax.hist(counts, bins='auto', log=True)
    ax.set_title(mech)
    ax.set_yscale("log")
    

for ax in axes[len(count_dict):]:
    ax.axis("off")

fig.supxlabel("Degree / Count")
fig.supylabel("Frequency")
fig.tight_layout()
plt.show()

# %%
