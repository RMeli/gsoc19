import MDAnalysis as mda
import pybel

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

import os

# PDBbind subsets 
datasets = ["refined", "other"]

# Dataset column names
col_names = ["n_modes", "n_flex"]

df = pd.DataFrame(columns=col_names, dtype=int)

no_conf=[] # Systems where no conformation was found in search space
failed=[] # Systems where a problem occurred

for dataset in datasets:
    for system in os.listdir(dataset):
        fname = os.path.join(f"{dataset}/{system}", "flex.pdb")
        
        try:
            u = mda.Universe(fname)

        except EOFError:
            print(f"No conformations found for system {system}.")
            no_conf.append(os.path.join(dataset, system))
            continue
        except Exception:
            print(f"Problem loading PDB file for system {system}")
            failed.append(os.path.join(dataset, system))
            continue

        n_modes = len(u.trajectory)
        n_flex = len(u.residues)

        data = pd.DataFrame([[n_modes, n_flex]], index=[system], columns=col_names)
        df = df.append(data)

print("\nDescription:")
print(df.astype(int).describe())

print(f"\nNo flexible residues printed ({len(failed)}):\n", *failed)
np.savetxt("analysis/noflex.dat", np.array(failed), "%s")

print(f"\nNo conformations found ({len(no_conf)}):\n", *no_conf)
np.savetxt("analysis/noconf.dat", np.array(no_conf), "%s")

# Number of rotatable bonds for systems where no conformation is found
n_rot_noconf = []
for nc in no_conf:
    path = os.path.join("../PDBbind18", nc)
    name = os.path.split(nc)[-1]
    ligand = os.path.join(path, name + "_ligand.mol2")

    for mol in pybel.readfile("mol2", ligand):
        n_rot_noconf.append(mol.OBMol.NumRotors())

print(f"\nNumber of rotable bond ({len(no_conf)}):\n", *n_rot_noconf)

# Plot histogram for number of modes
bins=np.arange(22) - 0.5 # Center-aligned bins
plt.hist(df["n_modes"].values, bins=bins, rwidth=0.9)
plt.xticks(range(21))
plt.xlim([-1,21])
plt.xlabel("Number of modes")
plt.savefig("analysis/plots/n_modes.pdf")
