import MDAnalysis as mda

import numpy as np
import pandas as pd

import os

# PDBbind subsets 
datasets = ["refined", "other"]

# Dataset column names
col_names = ["n_modes", "n_flex"]

df = pd.DataFrame(columns=col_names)

no_conf=[] # Systems where no conformation was found in search space
failed=[] # Systems where a problem occurred

for dataset in datasets:
    for system in os.listdir(dataset):
        fname = os.path.join(f"{dataset}/{system}", "flex.pdb")
        
        try:
            u = mda.Universe(fname)

        except EOFError:
            print(f"No conformations found for system {system}.")
            no_conf.append(system)
            continue
        except Exception:
            print(f"Problem loading PDB file for system {system}")
            failed.append(system)
            continue

        n_modes = len(u.trajectory)
        n_flex = len(u.residues)

        data = pd.DataFrame([[n_modes, n_flex]], index=[system], columns=col_names)
        df = df.append(data)

print("\nDescription:")
print(df.astype(float).describe())

print(f"\nFailed ({len(failed)}):\n", *failed)
np.savetxt("failed.dat", np.array(failed), "%s")

print(f"\nNo conformations found ({len(no_conf)}):\n", *no_conf)
np.savetxt("noconf.dat", np.array(no_conf), "%s")