import pandas as pd


file_path = '/Users/beatrizrodriguez/Desktop/conserved-quantities/mechanisms/cracmm/mech_cracmm3.def'
equations = []
current = ""
with open(file_path, "r") as f:
    for line in f:
        line = line.rstrip()
        if not line or line.lstrip().startswith("!"):
            continue
        current += " " + line.strip()
        if ";" in line:
            equations.append(current.strip())
            current = ""

clean_equations = []

for eq in equations:
    eq = eq.split("#")[0].strip()
    clean_equations.append(eq)


rows = []

for eq in clean_equations:
    rid, rest = eq.split(">", 1)
    rid = rid + ">"         

    reactants, products = rest.split("=", 1)

    rows.append({
        "reaction_id": rid.strip(),
        "reactants": reactants.strip(),
        "products": products.replace(";", "").strip()
    })

df = pd.DataFrame(rows)
df = df.reset_index()
df.loc[:, "index"] = df.loc[:, "index"]+1
df['index'] = 'R' + df['index'].astype(str)
df['reaction_id'] = df['reaction_id'].str.replace('<', '', regex=False)
df['reaction_id'] = df['reaction_id'].str.replace('>', '', regex=False)

orig_r_id = list(df.loc[:, "reaction_id"])

df = df.drop(columns=['reaction_id'])

d = {}

for r in range(df.shape[0]):
    row = df.iloc[r, :]
    idx = row[0]
    row[1] = row[1].replace("-", "+-")
    row[2] = row[2].replace("-", "+-")
    react_list = [] if pd.isna(row[1]) else row[1].strip().split("+")
    prod_list  = [] if pd.isna(row[2]) else row[2].strip().split("+")
    d[idx] = [react_list, prod_list]

species = set(
    elem
    for list_of_lists in d.values()
    for inner_list in list_of_lists
    for elem in inner_list
)

species = list({
    s.split("*")[-1].strip()
    for s in species
})

spc = {}
for i in range(len(species)):
    spc[species[i]] = i + 1

def parse_term(term):
    term = term.strip()
    if "*" in term:
        coeff_str, species = term.split("*", 1)
        coeff = float(coeff_str.replace(" ", ""))
    else:
        coeff = 1.0
        species = term
    return coeff, species.strip()

l = []
i = 0
for key, value in d.items():
    for reactant in value[0]:
        c1, r = parse_term(reactant)
        line = [orig_r_id[i], spc[r], key[1:], r, key, c1*-1]
        l.append(line)
    for product in value[1]:
        c2, p = parse_term(product)
        line2 = [orig_r_id[i], spc[p], key[1:], key, p, c2]
        l.append(line2)
    i+=1

df2 = pd.DataFrame(l, columns=["original r_id", "# species_index (starts from 1)","reaction_index (starts from 1)","from","to", " directed stoichiometric value"])

filename = "cracmm3_EdgeList.csv"
df2.to_csv(filename, index=False)
#%%

# No DELTA_C, DELTA_N, DELTA_SI
import pandas as pd


file_path = '/Users/beatrizrodriguez/Desktop/conserved-quantities/mechanisms/cracmm/mech_cracmm3.def'
equations = []
current = ""
with open(file_path, "r") as f:
    for line in f:
        line = line.rstrip()
        if not line or line.lstrip().startswith("!"):
            continue
        current += " " + line.strip()
        if ";" in line:
            equations.append(current.strip())
            current = ""

clean_equations = []

for eq in equations:
    eq = eq.split("#")[0].strip()
    clean_equations.append(eq)


rows = []

for eq in clean_equations:
    rid, rest = eq.split(">", 1)
    rid = rid + ">"         

    reactants, products = rest.split("=", 1)

    rows.append({
        "reaction_id": rid.strip(),
        "reactants": reactants.strip(),
        "products": products.replace(";", "").strip()
    })

df = pd.DataFrame(rows)
df = df.reset_index()
df.loc[:, "index"] = df.loc[:, "index"]+1
df['index'] = 'R' + df['index'].astype(str)
df['reaction_id'] = df['reaction_id'].str.replace('<', '', regex=False)
df['reaction_id'] = df['reaction_id'].str.replace('>', '', regex=False)


orig_r_id = list(df.loc[:, "reaction_id"])

df = df.drop(columns=['reaction_id'])

d = {}

for r in range(df.shape[0]):
    row = df.iloc[r, :]
    idx = row[0]
    row[1] = row[1].replace("-", "+-")
    row[2] = row[2].replace("-", "+-")
    react_list = [] if pd.isna(row[1]) else row[1].strip().split("+")
    prod_list  = [] if pd.isna(row[2]) else row[2].strip().split("+")
    d[idx] = [react_list, prod_list]

species = set(
    elem
    for list_of_lists in d.values()
    for inner_list in list_of_lists
    for elem in inner_list
)

species = list({
    s.split("*")[-1].strip()
    for s in species
})

eliminate = ["DELTA_C", "DELTA_N", "DELTA_SI"]

species_elim = [i for i in species if i not in eliminate]


spc = {}
for i in range(len(species_elim)):
    spc[species_elim[i]] = i + 1

def parse_term(term):
    term = term.strip()
    if "*" in term:
        coeff_str, species = term.split("*", 1)
        coeff = float(coeff_str.replace(" ", ""))
    else:
        coeff = 1.0
        species = term
    return coeff, species.strip()

l = []
i = 0
for key, value in d.items():
    for reactant in value[0]:
        c1, r = parse_term(reactant)
        if r in spc.keys():
            line = [orig_r_id[i], spc[r], key[1:], r, key, c1*-1]
            l.append(line)
    for product in value[1]:
        c2, p = parse_term(product)
        if p in spc.keys():
            line2 = [orig_r_id[i], spc[p], key[1:], key, p, c2]
            l.append(line2)
    i+=1

df2 = pd.DataFrame(l, columns=["original r_id", "# species_index (starts from 1)","reaction_index (starts from 1)","from","to", " directed stoichiometric value"])

filename = "cracmm3_EdgeList_eliminated.csv"
df2.to_csv(filename, index=False)
# %%
