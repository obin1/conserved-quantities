import pandas as pd
import re

from algorithm import create_sparse, create_coproduction, create_symbols, merge_coprod, s_linalg

#%%

edge_list_amore = pd.read_csv("../mechanisms/amore2_isop_135species/amore2_isop_135species_EdgeList.csv", comment="!")

# get unique values in first column (debugging the 135 species discrepancy)
unique_species = edge_list_amore.iloc[:, 0].unique()
print("number of unique species in edgelist: ", len(unique_species))

# double check, read in the reactions list and parse it
reactions = pd.read_csv("../mechanisms/amore2_isop_135species/amore2_isop_135species_Rnames.csv",header=None)
# parse reactions. split at + and =
species = []
for r in reactions.iloc[:,0]:
    parts = re.split(r'\s*\+\s*|\s*=\s*', r)
    for p in parts:
        name = p.strip()
        if not name:
            continue
        species.append(name)

print("number of unique species in rxn list: ", len(set(species)))