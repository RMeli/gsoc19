import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

sim = pd.read_csv("similarity.csv", index_col=0)

sim_mtx = sim.to_numpy()
assert sim_mtx.shape == (92, 92)
assert np.allclose(np.diag(sim_mtx), 1)
# assert np.allclose(sim_mtx, sim_mtx.T)

diff = sim.subtract(sim.transpose())


def nonzerocols(x):
    l = list(cols[x.values])
    return l if l else np.nan


cols = diff.columns
bt = diff.apply(lambda x: x > 0)  # Mask nonzero values
nosymm = bt.apply(nonzerocols, axis=1).dropna()  # Get column names for nonzero values
print(nosymm)

plt.figure()
sns.heatmap(sim)
plt.savefig("sim.png")

plt.figure()
sns.heatmap(sim.subtract(sim.transpose()))
plt.savefig("sim-T.png")

plt.figure()
sns.clustermap(sim)
plt.savefig("sim-cluster.png")

print(sim)
