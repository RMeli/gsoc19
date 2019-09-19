import MDAnalysis as mda

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

import os
import warnings

# PDBbind subsets 
datasets = ["refined", "other"]

# Dataset column names
col_names = ["n_modes", "n_flex"]

# Output directory
outdir = os.path.abspath("analysis")
plotdir = os.path.join(outdir, "plots")

df = pd.DataFrame(columns=col_names, dtype=int)

def load(fpath: str, print_warnings=False) -> mda.Universe:
    """
    Load file with MDAnalysis, suppressing warnings by default.

    MDAnslysis complains about MOL2 atom types and some metal atoms in PDB files.
    """

    if not print_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            return mda.Universe(fpath)

    return mda.Universe(fpath)

no_conf=[] # Systems where no conformation was found in search space
no_flex=[] # Systems where no flexible residues are found
failed=[] # Systems where a problem occurred

for dataset in datasets:
    for system in os.listdir(dataset):
        fname = os.path.join(f"{dataset}/{system}", "flex.pdb")
        
        try:
            u = load(fname)

        except EOFError: # Empty PDB file when no conformations found
            print(f"No conformations found for system {system}.")
            no_conf.append(os.path.join(dataset, system))
            continue

        except Exception: # PDB file with no residues or no PDB file
            if os.path.isfile(fname): # File exists
                print(f"No flexible residues found for system {system}.")
                no_flex.append(os.path.join(dataset, system))

                # Get number of modes from ligand poses
                dname = os.path.join(f"{dataset}/{system}", "dock.pdb")
                try: 
                    u = load(dname)

                    n_modes = len(u.trajectory)
                    n_flex = 0 # No flexible residues

                    data = pd.DataFrame([[n_modes, n_flex]], index=[system], columns=col_names)
                    df = df.append(data)

                except Exception:
                    print(f"Problem loading dock.pdb file for system {system}")
                    continue

            else: # File does not exist (docking failed before output)
                print(f"Problem loading flex.pdb file for system {system}")
                failed.append(os.path.join(dataset, system))
                continue

        n_modes = len(u.trajectory)
        n_flex = len(u.residues)

        data = pd.DataFrame([[n_modes, n_flex]], index=[system], columns=col_names)
        df = df.append(data)

if not os.path.exists(outdir):
    os.makedirs(outdir)

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

print("Description:")
print(df.astype(int).describe())

num_no_conf = len(no_conf)
if num_no_conf > 0:
    print(f"\nNo conformations found ({num_no_conf}):\n", *no_conf)
    np.savetxt(os.path.join(outdir, "noconf.dat"), np.array(no_conf), "%s")

num_no_flex = len(no_flex)
if num_no_flex > 0:
    print(f"\nNo flexible residues found ({num_no_flex}):\n", *no_flex)
    np.savetxt(os.path.join(outdir, "noflex.dat"), np.array(no_flex), "%s")

num_failed = len(failed)
if num_failed > 0:
    print(f"\nFailed({num_failed}):\n", *failed)
    np.savetxt(os.path.join(outdir, "failed.dat"), np.array(failed), "%s")

# Plot histogram for number of modes
plt.figure()
n = 20 # number of modes
bins=np.arange(n + 2) - 0.5 # Center-aligned bins
plt.hist(df["n_modes"].values, bins=bins, rwidth=0.9)
plt.xticks(range(n + 1))
plt.xlim([-1,n + 1])
plt.xlabel("Number of modes")
plt.ylabel("Number of systems")
plt.savefig(os.path.join(plotdir, "n_modes.pdf"))
plt.savefig(os.path.join(plotdir, "n_modes.png"))

# Plot histogram for number of flexible residues
plt.figure(figsize=(20,10))
n = df["n_flex"].max()
bins=np.arange(n + 2) - 0.5 # Center-aligned bins
plt.hist(df["n_flex"].values, bins=bins, rwidth=0.9)
plt.xticks(range(n + 1))
plt.xlim([-1,n + 1])
plt.xlabel("Number of flexible residues")
plt.ylabel("Number of systems")
plt.savefig(os.path.join(plotdir, "n_flexres.pdf"))
plt.savefig(os.path.join(plotdir, "n_flexres.png"))
