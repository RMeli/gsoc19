"""
Given a PDB index, load a protein-ligand complex from PDBBind, outlining residues
within FLEXDIST from the ligand.
"""

from pymol import cmd, stored

import os

def loadcomplex(system, dataset, flexdist=3, pdbbindpath="../PDBbind18"):

    # Clear everything
    cmd.reinitialize("everything")

    print(f"Loading {dataset}/{system}")

    ligand = os.path.join(pdbbindpath, dataset, system, f"{system}_ligand.mol2")
    receptor = os.path.join(pdbbindpath, dataset, system, f"{system}_protein.pdb")

    # Load ligand and receptor
    cmd.load(ligand, "ligand") # Selection name: ligand
    cmd.load(receptor, "receptor") # Selection name: receptor

    # Center and zoom to ligand
    cmd.center("ligand")
    cmd.zoom("ligand", 8)

    # Remove solvent
    cmd.remove("solvent")

    # Receptor rendering
    cmd.color("grey", "receptor")

    # Get residue index of atoms within FLEXDIST from the ligand
    noflex = ["ALA", "GLY", "PRO"]
    recsel = "receptor and not hydro " +  " ".join([f"and not resn {resn}" for resn in noflex])
    stored.list = []
    cmd.iterate(
        f"({recsel}) within {flexdist} of ligand", # Selection
        "stored.list.append((resn, resi, chain))" # Action
    )

    # Remove redundancies
    flexres = set(stored.list) # Set of flexible residues

    # Outline flexible residues
    for resn, resi, chain in flexres:
        if chain != "": # ???
            sel = f"receptor and (resi {resi} in chain {chain})"
            cmd.show("licorice", sel)
            cmd.color("red", sel)
            cmd.remove(f"hydro in ({sel})")

cmd.extend("loadcomplex", loadcomplex)