"""
Given a PDB index and a PDBbind DATASET ("refined" or "other"), load the following
structures:
- Original crystal structure (from PDBBINDPATH)
- All poses for the flexible residues (from flex.pdb in DOCKINGPATH)
- All poses for the ligand (from dock.pdb in DOCKINGPATH)
Then, show the following:
- Original crystal structure of the receptor (cartoon, gray)
- Original residues treated as flexible during flexible docking (licorice, gray)
- Ligand pose(s) of rank(s) IDX(s) (licorice, COLORS[IDX])
- Flexible residues pose(s) of rank(s) IDX(s) (licorice, COLORS[IDX])
"""

from pymol import cmd, stored

import os

num_modes = 20 # Number of docking modes

# Colors corresponding to docking rank
colors = [
    "chocolate", "skyblue", "limegreen", "warmpink", "limon", "violet",
    "brightorange", "sand", "lime", "deepteal", "hotpink", "yellowirange",
    "violepurple", "marine", "olive", "smudge", "deepsalmon", "splitpea",
    "lightteal", "slate",
]

def loadflexdock(system, dataset, idxs=["1"], flexdist="3", pdbbindpath="../PDBbind18", dockingpath=""):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idxs = idxs.split()
    idxs = [int(idx) for idx in idxs]
    if len(idxs) > 7:
        raise RuntimeError("Displaying more than 5 poses is not supported.")

    # Convert flexdist to float
    flexdist = float(flexdist)

    # Print informations
    print(f"Loading {dataset}/{system} (rank {idxs})")
    print(f"flexdist = {flexdist}")

    # Build paths
    ligandpath = os.path.join(dockingpath, dataset, system, f"dock.pdb") # Docked ligands
    flexrespath = os.path.join(dockingpath, dataset, system, f"flex.pdb") # Flexible residues
    cligandpath = os.path.join(pdbbindpath, dataset, system, f"{system}_ligand.mol2") # Crystal ligand
    receptorpath = os.path.join(pdbbindpath, dataset, system, f"{system}_protein.pdb") # Crystal receptor

    # Load ligand and receptor
    # When loading a multi-MODEL PDB file, PyMol appends 000X to the selection name
    # DOCKSEL and FLEXSEL dictionaries maps an index to the actual selection
    cmd.load(ligandpath, "ligand") # Selection name: ligand_0001, ligand_0002, ...
    cmd.load(cligandpath, "cligand") # Selection name: ligand_0001, ligand_0002, ...
    cmd.load(flexrespath, "flexres") # Selection name: flexres_0001, flexres_0002, ...
    cmd.load(receptorpath, "receptor") # Selection name: receptor

    # Map selection name to rank (0 is the crystal)
    docksel = {i : f"ligand_{i:04d}" for i in range(1,num_modes+1)}
    docksel[0] = "cligand" # Add index 0 for crystal
    flexsel = {i : f"flexres_{i:04d}" for i in range(1,num_modes+1)}
    flexsel[0] = f"(receptor and not hydro) within {flexdist} of {docksel[0]}" # Add index 0 for crystal
    
    # Hide everything and show only receptor and ligands
    cmd.hide("all")
    cmd.show("cartoon", "receptor")
    for idx in idxs:
        # Show ligand as ball and stick
        cmd.show("sticks", docksel[idx])
        cmd.show("spheres", docksel[idx])
        cmd.set("sphere_scale", 0.2, docksel[idx])
        
        # Show flexible residue as licorice
        cmd.show("licorice", flexsel[idx])
        cmd.set("stick_radius", 0.2, flexsel[idx])

        # Color C atoms of ligand and flexible residues
        cmd.color(colors[idx], docksel[idx] + " and name C*")
        cmd.color(colors[idx], flexsel[idx] + " and name C*")

    # Center and zoom to ligand
    if len(idxs) == 1: # Center on the single ligand
        cmd.center(docksel[idxs[0]])
        cmd.zoom(docksel[idxs[0]], 10)
    else: # Center on the crystal ligand if more than one ligand
        cmd.center(docksel[0])
        cmd.zoom(docksel[0], 8)

    # Remove solvent
    cmd.remove("solvent")

    # Receptor rendering
    cmd.color("grey", "receptor")

    # Get residue index of atoms within FLEXDIST from the crystal ligand
    # Exclude ALA, GLY, PRO 
    noflex = ["ALA", "GLY", "PRO"]
    recsel = "receptor and not hydro " +  " ".join([f"and not resn {resn}" for resn in noflex])
    stored.list = []
    cmd.iterate(
        f"(receptor and not hydro) within {flexdist} of {docksel[0]}", # Selection
        "stored.list.append((resn, resi, chain))" # Action
    )

    # Remove redundancies
    flexres = set(stored.list) # Set of flexible residues

    # Outline flexible residues
    for resn, resi, chain in flexres:
        if chain != "": # ???
            sel = f"receptor and (resi {resi} in chain {chain})"
            cmd.show("licorice", sel)
            cmd.color("grey", sel)
            cmd.remove(f"hydro in ({sel})")
            cmd.set("stick_radius", 0.2, sel)

    # Show metal atoms
    cmd.show("spheres", "metals")
    cmd.set("sphere_scale", 0.5, "metals")


cmd.extend("loadflexdock", loadflexdock)