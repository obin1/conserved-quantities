import pandas as pd

df = pd.read_csv("/Users/beatrizrodriguez/Desktop/conserved-quantities/mechanisms/cracmm/cracmm2_rxn_metadata.csv")

df = df.drop(columns=['reaction_family', 'reaction_phase', 'publication_string',
       'publication_doi', 'underlying_data_publication_string',
       'underlying_data_publication_doi', 'Notes'])

df = df.reset_index()
df[["index"]] = df[["index"]] + 1
df['index'] = 'R' + df['index'].astype(str)

d = {}

for r in range(df.shape[0]):
    row = df.iloc[r, :]
    idx = row[0]
    react_list = [] if pd.isna(row[2]) else row[2].strip().split("+")
    prod_list  = [] if pd.isna(row[3]) else row[3].strip().split("+")
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
        coeff = float(coeff_str)
    else:
        coeff = 1.0
        species = term
    return coeff, species

l = []
for key, value in d.items():
    for reactant in value[0]:
        c1, r = parse_term(reactant)
        line = [spc[r], key[1:], r, key, c1*-1]
        l.append(line)
    for product in value[1]:
        c2, p = parse_term(product)
        line2 = [spc[p], key[1:], key, p, c2]
        l.append(line2)

df2 = pd.DataFrame(l, columns=["# species_index (starts from 1)","reaction_index (starts from 1)","from","to", "directed stoichiometric value"
])

filename = "cracmm2_EdgeList.csv"
df2.to_csv(filename, index=False)