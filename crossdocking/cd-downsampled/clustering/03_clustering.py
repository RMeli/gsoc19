import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

sim = pd.read_csv("similarity.csv", index_col=0)

assert sim.to_numpy().shape == (92, 92)
#assert np.allclose(np.diag(sim.to_numpy()), 1)

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