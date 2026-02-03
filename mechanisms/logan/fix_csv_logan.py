import pandas as pd
import numpy as np

df = pd.read_csv("../mechanisms/logan/Logan1981_edgelist.csv")

df_list = []
r_id = []
spc_delete = ["aerosol","precip_loss_highaltitude","precip_loss_lowaltitude","M","N2","O2","hv","products_R25","H2O2_aerosol_loss"]

# if reactant multiply coefficient of stoich by -1
# take note of reaction IDs
for i in range(df.shape[0]):
    row = df.loc[i, :]
    if row[0].startswith("R"):
        if row[1] not in spc_delete:
            df_list.append(row)
            r_id.append(row[0])
    elif row[1].startswith("R"):
        if row[0] not in spc_delete:
            row[2] = row[2]*-1
            df_list.append(row)
            r_id.append(row[1])

# set of reaction IDs
r_id_clean = list(set(r_id))

# reaction ID column
# replace "O" with respective species
df2 = pd.DataFrame(df_list)
df2["R_ID"] = r_id
df2.loc[60, "to"] = "O3P"
df2.loc[139, "from"] = "O3P"
df2.loc[177, "to"] = "O1D"
df2.loc[189, "to"] = "O3P"
df2.loc[181, "to"] = "O3P"

# replace H2CO with CH2O
df2["from"] = df2["from"].str.replace("H2CO", "CH2O", regex=False)
df2["to"] = df2["to"].str.replace("H2CO", "CH2O", regex=False)

# list of species in "from"
spc = [i for i in list(df2.loc[:, "from"]) if not i.startswith("R")]
# list of species in "to"
spc2 = [i for i in df2.loc[:, "to"] if not i.startswith("R")]
# combine list of species, take the unique set
spc_mixed = spc + spc2
spc3 = np.unique(spc_mixed)
spc4 = np.setdiff1d(spc3,spc_delete)
# dictionary to map species to species index
spc_d = {}
# i = 1
for idx, spc in enumerate(spc4):
    spc_d[spc] = idx+1
    # i += 1

df2 = df2.reset_index()

# separate reactant df and product df
react = df2[df2["stoichiometry"] < 0].reset_index().drop(columns=["level_0"])
prod = df2[df2["stoichiometry"] > 0].reset_index().drop(columns=["level_0"])

# create list of species index for reactants
react_spc = []
for i in range(react.shape[0]):
    idx = react.loc[i, "from"]
    if idx in spc_d.keys():
        react_spc.append(spc_d[idx])
# create list of species index for products
prod_spc = []
for i in range(prod.shape[0]):
    idx = prod.loc[i, "to"]
    if idx in spc_d.keys():
        prod_spc.append(spc_d[idx])

# create species index column
react["spc_idx"] = react_spc
prod["spc_idx"] = prod_spc
# concat reactant and product dfs
df3 = pd.concat([react, prod])
# sort values to order by reaction ID again
sorted_df = df3.sort_values(by='index')
sorted_df["R_idx"] = sorted_df["R_ID"].str.replace("R", "", regex=False)

x = list(sorted_df["R_idx"])
ordered_x = []
for item in x:
    if item not in ordered_x:
        ordered_x.append(item)

r_dict = {}
for index, r in enumerate(ordered_x):
    r_dict[r] = index + 1

r_int = []
for r in x:
    r_int.append(r_dict[r])

sorted_df["R_idx"] = r_int

# reorder columns
new_order = ["index", "R_ID", "spc_idx","R_idx", "from", "to", "stoichiometry"]
df4 = sorted_df[new_order]

# rename columns
df4.rename(columns={'R_ID': 'original r_id', 'spc_idx': '# species_index (starts from 1)', "R_idx":"reaction_index (starts from 1)", "stoichiometry":" directed stoichiometric value"}, inplace=True)

df4 = df4.reset_index()
df4 = df4.drop(columns=["index", "level_0"])


filename = "log81_EdgeList.csv"
df4.to_csv(filename, index=False)