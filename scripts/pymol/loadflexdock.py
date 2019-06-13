from pymol import cmd, stored

import os

num_modes = 20 # Number of docking modes

def loadflexdock(system, dataset, idx=1, flexdist=3, pdbbindpath="../PDBbind18", dockingpath=""):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string idx to number
    idx = int(idx)
    flexdist = float(flexdist)

    print(f"Loading {dataset}/{system} (rank {idx})")
    print(f"flexdist = {flexdist}")

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
    flexsel[0] = f"receptor within {flexdist} of {docksel[idx]}" # Add index 0 for crystal
    
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
    
    # Remove hydrogen atoms
    #cmd.remove("hydro")

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