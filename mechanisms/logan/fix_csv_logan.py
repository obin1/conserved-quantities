import pandas as pd

df = pd.read_csv("../mechanisms/logan/Logan1981_edgelist.csv")

df_list = []
r_id = []

# if reactant multiply coefficient of stoich by -1
# take note of reaction IDs
for i in range(df.shape[0]):
    row = df.loc[i, :]
    if row[0].startswith("R"):
        df_list.append(row)
        r_id.append(row[0])
    elif row[1].startswith("R"):
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
spc = [i for i in df2.loc[:, "from"] if not i.startswith("R")]
# list of species in "to"
spc2 = [i for i in df2.loc[:, "to"] if not i.startswith("R")]
# combine list of species, take the unique set
spc3 = set(spc + spc2)
# dictionary to map species to species index
spc_d = {}
i = 1
for spc in spc3:
    spc_d[spc] = i
    i += 1

df2 = df2.reset_index()

# separate reactant df and product df
react = df2[df2["stoichiometry"] < 0].reset_index().drop(columns=["level_0"])
prod = df2[df2["stoichiometry"] > 0].reset_index().drop(columns=["level_0"])

# create list of species index for reactants
react_spc = []
for i in range(react.shape[0]):
    idx = react.loc[i, "from"]
    react_spc.append(spc_d[idx])
# create list of species index for products
prod_spc = []
for i in range(prod.shape[0]):
    idx = prod.loc[i, "to"]
    prod_spc.append(spc_d[idx])

# create species index column
react["spc_idx"] = react_spc
prod["spc_idx"] = prod_spc
# concat reactant and product dfs
df3 = pd.concat([react, prod])
# sort values to order by reaction ID again
sorted_df = df3.sort_values(by='index')
sorted_df["R_idx"] = sorted_df["R_ID"].str.replace("R", "", regex=False)

# reorder columns
new_order = ["index", "R_ID", "spc_idx","R_idx", "from", "to", "stoichiometry"]
df4 = sorted_df[new_order]

# rename columns
df4.rename(columns={'R_ID': 'original r_id', 'spc_idx': '# species_index (starts from 1)', "R_idx":"reaction_index (starts from 1)", "stoichiometry":" directed stoichiometric value"}, inplace=True)
# remove unwanted species
df4 = df4[~((df4["from"] == "M") | (df4["to"] == "M"))]
df4 = df4[~((df4["from"] == "N2") | (df4["to"] == "N2"))]
df4 = df4[~((df4["from"] == "O2") | (df4["to"] == "O2"))]
df4 = df4[~((df4["from"] == "hv") | (df4["to"] == "hv"))]



df4 = df4.reset_index()
df4 = df4.drop(columns=["index", "level_0"])

# remove chlorine and wet deposition products
df4 = df4.drop(183)
df4 = df4.drop(180)
df4 = df4.drop(173)
df4 = df4.drop(86)
df4 = df4.reset_index()
df4 = df4.drop(columns=["index"])


filename = "log81_EdgeList.csv"
df4.to_csv(filename, index=False)