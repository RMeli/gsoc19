"""
Create training and test files with annotations.

Notes
-----
df.rank is a method of pd.DataFrame; the "rank" column need to be accessed with df["rank"]
"""

import os

import pandas as pd
import tqdm

import argparse

from sklearn.model_selection import GroupKFold

root = "wierbowski_cd"

def ligrmsd(row):
    if row.rmsd < 2.0:
        return 1  # Good pose
    else:
        return 0  # Bad pose

parser = argparse.ArgumentParser(description="Cross-validation folds.")
parser.add_argument(
    "-c",
    "--cluster",
    action="store_true",
    help="Output cluster index in affinity column.",
)
args = parser.parse_args()

prefix = "cluster" if args.cluster else "nc"

def ligname(row):
    fprefix = f"{row.protein}_PRO_{row.ligand}_LIG_aligned"
    fsuffix = f"_default_ensemble_none_flexdistNone_p{row['rank']}.gninatypes"

    for v in ["", "_v2"]:
        fname = os.path.join(root, row.pocket, fprefix + v + fsuffix)
        if os.path.isfile(fname):
            break

    return fname


def recname(row):
    return os.path.join(root, row.pocket, f"{row.protein}_PRO.gninatypes")

df = pd.read_csv("wierbowski_cd/rmsds.csv", index_col=False)

# Create CV clusters based on ProBiS
clusters = pd.read_csv("../cd-downsampled/clustering/clusters.csv", index_col=0)
def pocket_to_cluster(row):
    return clusters.loc[row.pocket]

df["group"] = df.apply(pocket_to_cluster, axis=1)

# Annotate DataFrame
df["annotation"] = df.apply(ligrmsd, axis=1)

# Detect empty groups (only decoys after annotation)
empty = []
for group, dfgroup in df.groupby(by="group"):
    if all(dfgroup["annotation"] == 0):
        print(f"Group {group} only contains decoys. Removing...")
        empty.append(group)

# Remove empty groups
for e in empty:
    df = df[df["group"] != e]

df.to_csv(f"analysis/{prefix}_annotated.csv", index=False, float_format="%.5f")

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
                line = f"{row.annotation} {row.group} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            else:
                line = f"{row.annotation} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"

            trout.write(line)

    with open(f"files/{prefix}test{fold}.types", "w") as teout:
        for te in tqdm.tqdm(test_idx, leave=False, desc=f"Test {fold}"):
            row = df.iloc[te]

            if args.cluster:
                line = f"{row.annotation} {row.group} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            else:
                line = f"{row.annotation} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"

            teout.write(line)
