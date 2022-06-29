import csv
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
plt.style.use("seaborn")


# question 2e):
rows = []
with open("gene_disease_opt.csv", 'r') as file:
    data = csv.reader(file)
    header = next(data)
    for row in data:
        rows.append(row)
# print(header)
# print(rows)

# question 2f):
# pd_data = pd.read_csv("gene_disease_opt.csv", sep=' ', on_bad_lines='skip')
pd_data = pd.read_csv("gene_disease_opt.csv", sep="\t", on_bad_lines='skip')
pd_data.sort_values(["gene_family"], axis=0, ascending=[False], inplace=True)
pd_data.drop(pd_data.index[pd_data['score'] < 0.4], inplace=True)
diseases = []
for el in pd_data["disease.name"]:
    if el not in diseases:
        diseases.append(el)
        if len(diseases) == 5:
            break
print(diseases)

result = pd.DataFrame(data=diseases, columns=["diseases"])
print(result)
result.to_csv("../output/selected_associations.csv")


# question 2g):
associations = 0
for el in pd_data["publicationYear"]:
    if el > 2000:
        associations += 1

# Number of publications after 2020
print(associations)

count = 0
for el in pd_data["variantFunctionalConsequence.label"]:
    if el == "missense_variant":
        count += 1
# This prints the amount of missense variants out of the publications after 2020
print(count)


# question 2h):
fig, axes = plt.subplots(3, 1)
pd_data_1 = pd.read_csv("gene_disease_opt.csv", sep="\t", on_bad_lines='skip')
axes[0].hist(pd_data_1["score"][pd_data_1["gene_family"]=="OTHER"], bins=20, density=True, color="C0", alpha=0.5)
axes[1].hist(pd_data_1["score"][pd_data_1["gene_family"]=="SLC"], bins=20, density=True, color="C1", alpha=0.5)
axes[2].hist(pd_data_1["score"][pd_data_1["gene_family"]=="GPCR"], bins=20, density=True, color="C5", alpha=0.5)
plt.show()

fig.savefig('../graphics/score_dst.pdf')

# question 2i):
def transformed_data(data):
    # df = pd.read_csv(data)
    new_data = data.iloc[:, [0, 6]]
    return new_data.sort_values(["variantId"])

# As I am working in python and an index matrix is more inherent to R - I thought
# to return a pandas dataframe as that too also indexes the entries of the rows.
# The below is just to test the function is working
data_2 = pd_data_1.copy()
print(transformed_data(data_2))

# question 2j):
df = pd.read_csv("gene_disease_opt.csv", sep="\t", on_bad_lines='skip')
G = nx.from_pandas_edgelist(df, source=df["gene_family"], target=df["disease.id"], edge_attr=df["score"])
widths = np.array([w for *_, w in G.edges.data('weight')])

pos = nx.spring_layout(G, seed=7)

# nodes
nx.draw_networkx_nodes(G, pos, node_size=700)

# edges
nx.draw_networkx_edges(G, pos, width=widths*10)

# labels
nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.show()


