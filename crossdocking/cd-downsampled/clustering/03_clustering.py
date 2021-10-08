import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

import networkx as nx

# Z score threshold
ZscoreT = 3.5

zscore = pd.read_csv("zscore.csv", index_col=0)
zscore = zscore.sort_index()  # Sort index
zscore = zscore[sorted(zscore.columns)]

# Symmtetize
# zscore = zscore.add(zscore.transpose()) / 2

plt.figure()
sns.heatmap(zscore, vmin=0, vmax=5)
plt.savefig("zscore.png")

sim_mtx = zscore.to_numpy(copy=True)
sim_mtx[sim_mtx < ZscoreT] = 0
sim_mtx[sim_mtx >= ZscoreT] = 1

print(sim_mtx)

assert sim_mtx.shape == (92, 92)
assert np.allclose(np.diag(sim_mtx), 1)
# assert np.allclose(sim_mtx, sim_mtx.T)

diff = pd.DataFrame(
    data=sim_mtx - sim_mtx.T, columns=zscore.columns, index=zscore.index
)


def nonzerocols(x):
    l = list(cols[x.values])
    return l if l else np.nan


cols = diff.columns
bt = diff.apply(lambda x: x > 0)  # Mask nonzero values
nosymm = bt.apply(nonzerocols, axis=1).dropna()  # Get column names for nonzero
print(nosymm)

# FIXME: Is this necessary given the graph approach below?
# FIXME: Using undeirected graphs, should automatically get rid of this problem...
sim = pd.DataFrame(sim_mtx, index=zscore.index, columns=zscore.columns)
for name, val in nosymm.iteritems():
    for v in val:
        zscore1 = zscore.loc[name, v]
        zscore2 = zscore.loc[v, name]
        print(
            f"zscore({name:5}|{v:5})={zscore1:.2f} <---> zscore({v:5}|{name:5})={zscore2:.2f}\t{(zscore1 + zscore2) / 2:.2f}"
        )

        # Given the high z score threshold, manually set the systems as similar
        sim.loc[name, v] = 1.0
        sim.loc[v, name] = 1.0

sim_mtx = sim.to_numpy()
assert np.allclose(np.diag(sim_mtx), 1)
assert np.allclose(sim_mtx, sim_mtx.T)

plt.figure()
sns.heatmap(sim)
plt.savefig("sim.png")

plt.figure()
cm = sns.clustermap(sim)
plt.savefig("sim-cluster.png")

assert np.allclose(cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind)

# Build graph using similarity matrix as adjacency matrix
G = nx.Graph()
for i, p in enumerate(sim.index):
    G.add_node(i, pocket=p)

for i, p1 in enumerate(sim.index):
    for j, p2 in enumerate(sim.columns):
        if sim.loc[p1, p2] > 0.5:  # Similar pockets
            G.add_edge(i, j)

# Check that the adjacency matrix encodes the similarity
A = nx.to_numpy_array(G)
assert np.allclose(A, sim.to_numpy())

# Get different indices for different clusters
clusters = -np.ones(92)
for cc_idx, cc in enumerate(nx.connected_components(G)):
    for node_idx in cc:
        # Assign same cc_idx (cluster index) to nodes in the same connected component
        clusters[node_idx] = cc_idx

assert (clusters > -1).all()

for i, c in enumerate(clusters):
    sim.iloc[i, :] *= c + 1

plt.figure()
cm = sns.heatmap(
    sim.iloc[cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind],
    cmap="cividis",
    mask=sim.iloc[
        cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind
    ].to_numpy()
    < 0.5,
)
plt.savefig("sim-cluster-clustered.png")
