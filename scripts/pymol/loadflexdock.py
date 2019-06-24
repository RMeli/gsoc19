"""
Given a PDB index and a PDBbind DATASET ("refined" or "other"), load the following
structures:
- Original crystal structure (from PDBBINDPATH)
- Reconstructed crystal structure with selected flexible residues (from DOCKINGPATH)
- Selected pose for the ligand (from DOCKINGPATH)
Then, show the following:
- Reconstructed structure of the receptor (cartoon, gray)
- Original residues treated as flexible during flexible docking (licorice, purple)
- Flexible residues pose of rank IDX (licorice, real)
- Flexible residues atoms in reconstructed receptor, including backbone (sphere, yellow)
- Ligand pose of rank IDX (licorice, marine)

The  purpose of this script is to asses the successful reconstruction of the receptor
by substitution of the flexible residues within the original crystal structure.
"""

from pymol import cmd, stored

import os


def loadflexdock(
    system,
    dataset,
    idx="1",
    flexdist="3",
    pdbbindpath="../../PDBbind18",
    dockingpath="",
):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idx = int(idx)

    # Convert flexdist to float
    flexdist = float(flexdist)

    # Print informations
    print(f"Loading {dataset}/{system} (rank {idx})")
    print(f"flexdist = {flexdist}")

    # Build paths
    ligname = f"{system}_ligand-{idx}.pdb"
    ligandpath = os.path.join(dockingpath, dataset, system, ligname)  # Docked ligand
    flexname = f"{system}_flex-{idx}.pdb"
    flexpath = os.path.join(dockingpath, dataset, system, flexname)  # Flexible residues
    receptorname = f"{system}_protein-{idx}.pdb"
    receptorpath = os.path.join(dockingpath, dataset, system, receptorname)
    crystalname = f"{system}_protein.pdb"
    crystalpath = os.path.join(
        pdbbindpath, dataset, system, crystalname
    )  # Crystal receptor

    # Load ligand and receptor
    cmd.load(ligandpath, "ligand")  # Selection name: ligand
    cmd.load(flexpath, "flex")  # Selection name: flex
    cmd.load(receptorpath, "receptor")  # Selection name: receptor
    cmd.load(crystalpath, "crystal")  # Selection name: receptor

    # Hide everything
    cmd.hide("all")

    # Show receptor
    cmd.show("cartoon", "receptor")
    cmd.color("grey", "receptor")

    # Show docked ligand
    cmd.show("sticks", "ligand")
    cmd.show("spheres", "ligand")
    cmd.set("sphere_scale", 0.2, "ligand")
    cmd.color("marine", "ligand and name C*")

    # Show flexible residues
    cmd.show("licorice", "flex")
    cmd.set("stick_radius", 0.2, "flex")
    cmd.color("teal", "flex and name C*")

    # Center and zoom to ligand
    cmd.center("ligand")
    cmd.zoom("ligand", 10)

    # Remove solvent
    cmd.remove("solvent")

    # Get residues of receptor close to flexible residues
    stored.list = []
    cmd.iterate(
        "receptor within 0.1 of flex",  # Selection
        "stored.list.append((resn, resi, chain))",  # Action
    )

    # Remove redundancies
    flexres = set(stored.list)  # Set of flexible residues

    # Outline flexible residues of the receptor
    for _, resi, chain in flexres:
        if chain != "":  # ???
            sel = f"receptor and (resi {resi} in chain {chain})"
            cmd.show("sphere", sel)
            cmd.set("sphere_scale", 0.2, sel)
            cmd.color("yellow", sel)
            cmd.remove(f"hydro in ({sel})")

    # Outline flexible residues of the crystal
    for _, resi, chain in flexres:
        if chain != "":  # ???
            sel = f"crystal and (resi {resi} in chain {chain})"
            cmd.show("licorice", sel)
            cmd.set("stick_radius", 0.15, sel)
            cmd.color("deeppurple", sel)
            cmd.remove(f"hydro in ({sel})")

    # Show metal atoms
    cmd.show("spheres", "metals")
    cmd.set("sphere_scale", 0.5, "metals")


cmd.extend("loadflexdock", loadflexdock)
