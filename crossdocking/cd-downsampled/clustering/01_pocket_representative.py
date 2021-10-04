"""
Similarity between pockets uses a representative structure.

Here, the representative sturcure is obtained using the cognate receptor of the
largest ligand. The largest ligand is selected by molecular weight.
"""
import os
import re
import pickle

from rdkit import Chem
from rdkit.Chem import Descriptors

dataroot = "../carlos_cd"

# Check that all pockets have been assigned to a cluster
pockets = [d for d in os.listdir(dataroot) if os.path.isdir(os.path.join(dataroot, d))]
# print(pockets)
assert len(pockets) == 92  # Check there are exactly 92 pockets

# List of systems
# (pocket_dir, ligname)
systems = []


for pocket in pockets:
    d = os.path.join(dataroot, pocket, "PDB_Structures")
    ligands = [
        flig for flig in os.listdir(d) if re.search("^.{4}_LIG_aligned.*\.sdf$", flig)
    ]
    # print(pocket, ligands)

    wmax = -1
    lwmax = ""
    for flig in ligands:
        supplier = Chem.SDMolSupplier(
            os.path.join(dataroot, pocket, "PDB_Structures", flig)
        )
        mol = next(supplier)

        # Select largest ligand in the pocket
        mw = Descriptors.MolWt(mol)
        # print(flig, mw)
        if mw > wmax:
            wmax = mw
            lwmax = flig

    systems.append((pocket, lwmax[:4]))

# print(systems)

with open("representatives.list", "w") as fout:
    for system in systems:
        fout.write(f"{system[0]}\t{system[1]}\n")
    # pickle.dump(systems, fout)
