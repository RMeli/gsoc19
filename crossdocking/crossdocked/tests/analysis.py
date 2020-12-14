import prody

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import numpy as np

from tqdm import tqdm

from matplotlib import pyplot as plt
import seaborn as sns

import os
import sys
import gzip
import re

# Turn of logging
prody.confProDy(verbosity="none")
rdkit.RDLogger.DisableLog('rdApp.error')

# Root folder for analysis
# root = sys.argv[1]

TIMINGS = {}
N_POSES = {}
N_ROTATABLE_BONDS = {}
N_FLEXIBLE_RESIDUES = {}

for root in ["docking-d3.0", "docking-d3.5"]:
    pockets = os.listdir(os.path.join(root, "docking"))

    timings = []
    n_poses = []
    n_rotatable_bonds = []
    n_flexible_residues = []

    for pocket in tqdm(pockets):
        # Path to pocket
        ppath = os.path.join(root, "docking", pocket)

        runs = [
            f.replace("-gnina.log", "") for f in os.listdir(ppath) if f.endswith(".log")
        ]

        for run in runs:
            log = f"{run}-gnina.log"
            lig = f"flexlig-{run}.sdf.gz"
            flx = f"flexrec-{run}.pdb.gz"

            with open(os.path.join(ppath, log), "r") as flog:
                content = flog.read()

                if re.search(
                    "WARNING: Could not find any conformations" +
                    "completely within the search space.",
                    content,
                ):
                    continue
                else:
                    target = "Loop time"

                    try:
                        time = int(
                                re.search(f"{target} \d+", content)
                                .group(0)
                                .replace(target, "")
                                .strip()
                            )
                    except AttributeError: # Problem loading the system
                        continue
                    else:
                        sdf = Chem.ForwardSDMolSupplier(gzip.open(os.path.join(ppath, lig), "r"))

                        poses = 0
                        for mol in sdf:
                            if mol is None:
                                # Some molecules can't be correctly parsed
                                continue

                            if poses == 0:  # Count same molecule only once
                                n_rotatable_bonds.append(
                                    rdMolDescriptors.CalcNumRotatableBonds(mol)
                                )

                                timings.append(time)

                                model = prody.parsePDB(os.path.join(ppath, flx), model=1)
                                try:
                                    n_flexible_residues.append(model.numResidues())
                                except AttributeError:
                                    n_flexible_residues.append(0)

                            poses += 1

                        n_poses.append(poses)

                        

    TIMINGS[root] = timings
    N_POSES[root] = n_poses
    N_ROTATABLE_BONDS[root] = n_rotatable_bonds
    N_FLEXIBLE_RESIDUES[root] = n_flexible_residues


def distplot(data, name, xlabel, bins=True):
    plt.figure()

    m, M = np.Inf, -np.Inf
    for _, value in data.items():
        if min(value) < m:
            m = min(value)
        if max(value) > M:
            M = max(value)

    b = np.arange(m - 0.5, M + 1, 1)

    for key, value in data.items():
        sns.distplot(np.array(value), kde=False, bins=b, label=key)

    plt.xticks(b[::2] + 0.5)

    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(f"{name}.pdf")
    plt.savefig(f"{name}.png")


def timeplot(timings, xvalues, xlabel):
    plt.figure()

    for key, times in timings.items():
        sns.scatterplot(x=xvalues[key], y=np.array(times) / 60.0, label=key)

    plt.xlabel(xlabel)
    plt.ylabel("time (min)")

    plt.yscale("log")
    plt.legend()

    plt.savefig(f"timings-{xlabel}.pdf")
    plt.savefig(f"timings-{xlabel}.png")


timeplot(TIMINGS, N_ROTATABLE_BONDS, "n_flex_residues")
timeplot(TIMINGS, N_FLEXIBLE_RESIDUES, "n_rot_bonds")
distplot(N_POSES, "n_poses", "Poses")
distplot(N_ROTATABLE_BONDS, "n_rotatable_bonds", "Rotatable bonds")
distplot(N_FLEXIBLE_RESIDUES, "n_flexible_residues", "Flexible residues")
