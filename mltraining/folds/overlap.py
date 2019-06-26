import pandas as pd
import os
import re

csvfile = "probis.csv"

df = pd.read_csv(csvfile)

pdbs = df["pdb"].to_numpy()

pdbbind = "../../PDBbind18"

refinedpath = os.path.join(pdbbind, "refined")
refined = [
    pdb for pdb in os.listdir(refinedpath)
    if re.match("^....$", pdb) and os.path.isdir(os.path.join(refinedpath, pdb))
]

found = 0
for pdb in pdbs:
    if pdb in refined:
        found += 1
print(f"Refined: found {found} out of {len(refined)}")


otherpath = os.path.join(pdbbind, "other")
other = [
    pdb for pdb in os.listdir(otherpath)
    if re.match("^....$", pdb) and os.path.isdir(os.path.join(otherpath, pdb))
]

found = 0
for pdb in pdbs:
    if pdb in other:
        found += 1
print(f"Other: found {found} out of {len(other)}")

