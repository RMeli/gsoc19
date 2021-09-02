"""

Notes
-----
df.rank is a method of pd.DataFrame; the "rank" column need to be accessed with df["rank"]
"""

import os

import pandas as pd

from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import OrdinalEncoder

root="carlos-cd"

def ligname(row):
    fprefix = f"{row.protein}_PRO_{row.ligand}_aligned"
    fsuffix = f"_default_ensemble_none_flexdist3.5_p{row['rank']}.sdf.gz"

    for v in ["", "_v2"]:
        fname = os.path.join(root, row.pocket, "PDB_Structures", fprefix + v + fsuffix)
        if os.path.isfile(fname):
            break

    return fname

def recname(row):
    fprefix = f"{row.protein}_PRO_{row.ligand}_aligned"
    fsuffix = f"_default_ensemble_none_flexdist3.5_full_p{row['rank']}.pdb.gz"

    for v in ["", "_v2"]:
        fname = os.path.join(root, row.pocket, "PDB_Structures", fprefix + v + fsuffix)
        if os.path.isfile(fname):
            break

    return fname

def ligrmsd(row):
    if row.rmsd < 2.0:
        return 1 # Good pose
    else:
        return 0 # Bad pose


df = pd.read_csv("analysis/rmsd_clean.csv", index_col=False)

# Assign different group to different pockets
# TODO: Assesss pocket similarity?
encoder = OrdinalEncoder()
df["group"] = encoder.fit_transform(df["pocket"].to_numpy().reshape(-1, 1))
df["group"] = df["group"].astype(int)

# Check that the number of unique pockets matches the number of groups
assert len(df["pocket"].unique()) == df["group"].max() + 1

# Split dataset in three folds (cross-validation)
# Randomly shuffle the data and keep the same pocket in the same fold
cv = GroupShuffleSplit(n_splits=3, random_state=42)
for fold, (train_idx, test_idx) in enumerate(
    cv.split(df.drop(columns=["group"]).to_numpy(), groups=df["group"].to_numpy())
):
    with open(f"files/train_{fold}.types", "w") as trout:
        for tr in train_idx:
            row =  df.iloc[tr]

            a = ligrmsd(row)
            line=f"{a} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            trout.write(line)

    with open(f"files/test_{fold}.types", "w") as teout:
        for te in test_idx:
            row =  df.iloc[te]

            a = ligrmsd(row)
            line=f"{a} {recname(row)} {ligname(row)} # {row.rmsd:.4f} {row.score:.4f}\n"
            teout.write(line)

