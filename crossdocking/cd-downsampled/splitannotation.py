"""
Create training and test files with split annotation for ligand pose and receptor pose.

Notes
-----
df.rank is a method of pd.DataFrame; the "rank" column need to be accessed with df["rank"]
"""

import os

import pandas as pd
import tqdm

import argparse

from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import OrdinalEncoder

root = "carlos_cd"

# Path to CSV file with RMSDs
csv_file = "analysis/rmsd2_clean.csv"


def ligrmsd(row):
    if row.rmsd <= 2.0:
        return 1  # Good pose
    else:
        return 0  # Bad pose


def recrmsd(row, col, t):
    if row[col] <= t:
        return 1  # Good pose
    else:
        return 0  # Bad pose


# Define all possible annotations
annotations = {
    "flex1": lambda r: recrmsd(r, "flexrmsd", 1.0),
    "flex2": lambda r: recrmsd(r, "flexrmsd", 2.0),
    "max2": lambda r: recrmsd(r, "fmaxrmsd", 2.0),
}

parser = argparse.ArgumentParser(description="Cross-validation folds.")
parser.add_argument(
    "-c",
    "--cluster",
    action="store_true",
    help="Output cluster index in affinity column.",
)
parser.add_argument(
    "--receptor",
    default=None,
    help="Receptor annotation",
    choices=annotations.keys(),
)
args = parser.parse_args()

prefix = "cluster" if args.cluster else "nc"


def ligname(row):
    fprefix = f"{row.protein}_PRO_{row.ligand}_LIG_aligned"
    fsuffix = f"_default_ensemble_none_flexdist3.5_p{row['rank']}.gninatypes"

    for v in ["", "_v2"]:
        fname = os.path.join(root, row.pocket, "PDB_Structures", fprefix + v + fsuffix)
        if os.path.isfile(fname):
            break

    return fname


def recname(row):
    fprefix = f"{row.protein}_PRO_{row.ligand}_LIG_aligned"
    fsuffix = f"_default_ensemble_none_flexdist3.5_full_p{row['rank']}.gninatypes"

    for v in ["", "_v2"]:
        fname = os.path.join(root, row.pocket, "PDB_Structures", fprefix + v + fsuffix)
        if os.path.isfile(fname):
            break

    return fname


df = pd.read_csv(csv_file, index_col=False)

# Create CV clusters based on ProBiS
clusters = pd.read_csv("clustering/clusters.csv", index_col=0)
# print(clusters)
def pocket_to_cluster(row):
    return clusters.loc[row.pocket]


df["group"] = df.apply(pocket_to_cluster, axis=1)

# Ligand annotation
df["lig_annotation"] = df.apply(ligrmsd, axis=1)

# Receptor annotation
df["rec_annotation"] = df.apply(annotations[args.receptor], axis=1)

prefix = "SPLIT" + args.receptor + prefix

# Detect empty groups (only decoys after annotation)
empty = []
for group, dfgroup in df.groupby(by="group"):
    if all(dfgroup["lig_annotation"] == 0):
        print(f"Group {group} only contains decoys. Removing...")
        empty.append(group)

# Remove empty groups
for e in empty:
    df = df[df["group"] != e]

df.to_csv(f"analysis/{prefix}_annotated.csv", index=False, float_format="%.5f")

# print(df)
# print(df["group"].unique())

# Split dataset in three folds (cross-validation)
# Randomly shuffle the data and keep the same pocket in the same fold
cv = GroupKFold(n_splits=3)
for fold, (train_idx, test_idx) in enumerate(
    cv.split(df.drop(columns=["group"]).to_numpy(), groups=df["group"].to_numpy())
):
    with open(f"files/{prefix}train{fold}.types", "w") as trout:
        for tr in tqdm.tqdm(train_idx, leave=False, desc=f"Train {fold}"):
            row = df.iloc[tr]

            if args.cluster:
                line = f"{row.lig_annotation} {row.rec_annotation} {row.group} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            else:
                line = f"{row.lig_annotation} {row.rec_annotation} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"

            trout.write(line)

    with open(f"files/{prefix}test{fold}.types", "w") as teout:
        for te in tqdm.tqdm(test_idx, leave=False, desc=f"Test {fold}"):
            row = df.iloc[te]

            if args.cluster:
                line = f"{row.lig_annotation} {row.rec_annotation} {row.group} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            else:
                line = f"{row.lig_annotation} {row.rec_annotation} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"

            teout.write(line)
