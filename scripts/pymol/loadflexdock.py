from pymol import cmd, stored

import os

num_modes = 20 # Number of docking modes

colors = [("blue", "skyblue"), ("red", "ruby"), ("green", "lime")]

def loadflexdock(system, dataset, idxs=["1"], flexdist="3", pdbbindpath="../PDBbind18", dockingpath=""):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idxs = idxs.split()
    idxs = [int(idx) for idx in idxs]
    if len(idxs) > 3:
        raise RuntimeError("Displaying more than 3 poses is not supported.")

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
    for i, idx in enumerate(idxs):
        cmd.show("licorice", docksel[idx])
        cmd.color(colors[i][0], docksel[idx])
        cmd.show("licorice", flexsel[idx])
        cmd.color(colors[i][1], flexsel[idx])

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

cmd.extend("loadflexdock", loadflexdock)