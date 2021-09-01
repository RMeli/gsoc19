"""

Notes
-----
df.rank is a method of pd.DataFrame; the "rank" column need to be accessed with df["rank"]
"""

import os

import pandas as pd

df = pd.read_csv("analysis/rmsd_clean.csv")

root = "carlos_cd"

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

for idx, row in df.iterrows():
    annotation = ligrmsd(row)
    lig = ligname(row)
    rec = recname(row)

    print(annotation, lig, rec)
    #break