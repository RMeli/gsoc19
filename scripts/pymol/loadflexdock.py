from pymol import cmd, stored

import os

num_modes = 20 # Number of docking modes

docksel = {f"{i}" : f"ligand_{i:04d}" for i in range(1,num_modes+1)}
flexsel = {f"{i}" : f"flexres_{i:04d}" for i in range(1,num_modes+1)}

def loadflexdock(system, dataset, idx=1, flexdist=3, pdbbindpath="../PDBbind18", dockingpath=""):

    print(f"Loading {dataset}/{system}")

    ligandpath = os.path.join(dockingpath, dataset, system, f"dock.pdb")
    flexrespath = os.path.join(dockingpath, dataset, system, f"flex.pdb")
    receptorpath = os.path.join(pdbbindpath, dataset, system, f"{system}_protein.pdb")

    # Load ligand and receptor
    # When loading a multi-MODEL PDB file, PyMol appends 000X to the selection name
    # DOCKSEL and FLEXSEL dictionaries maps an index to the actual selection
    cmd.load(ligandpath, "ligand") # Selection name: ligand_0001, ligand_0002, ...
    cmd.load(flexrespath, "flexres") # Selection name: flexres_0001, flexres_0002, ...
    cmd.load(receptorpath, "receptor") # Selection name: receptor
    
    # Hide everything and show only receptor and single ligand
    cmd.hide("all")
    cmd.show("cartoon", "receptor")
    cmd.show("licorice", docksel[idx])
    cmd.show("licorice", flexsel[idx])

    # Center and zoom to ligand
    cmd.center(docksel[idx])
    cmd.zoom(docksel[idx], 8)

    # Remove solvent
    cmd.remove("solvent")

    # Receptor rendering
    cmd.color("grey", "receptor")

    # Get residue index of atoms within FLEXDIST from the ligand
    stored.list = []
    cmd.iterate(
        f"all within {flexdist} of {docksel[idx]}", # Selection
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

cmd.extend("loadflexdock", loadflexdock)