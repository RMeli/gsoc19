from pymol import cmd, stored

import os, sys

sys.path.insert(1, '../../scripts/python/')

import flexrmsd as frmsd

def flexrmsd(
    system,
    dataset,
    idx="1",
    flexdist="3",
    pdbbindpath="../../PDBbind18",
):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idx = int(idx)

    # Convert flexdist to float
    flexdist = float(flexdist)

    # Print informations
    print(f"Loading {dataset}/{system} (rank {idx})")

    # Build paths
    ligname = f"{system}_ligand-{idx}.pdb"
    ligandpath = os.path.join(dataset, system, ligname)  # Docked ligand
    flexname = f"{system}_flex-{idx}.pdb"
    flexpath = os.path.join(dataset, system, flexname)  # Flexible residues
    receptorname = f"{system}_protein-{idx}.pdb"
    receptorpath = os.path.join(dataset, system, receptorname)
    crystalname = f"{system}_protein.pdb"
    crystalpath = os.path.join(
        pdbbindpath, dataset, system, crystalname
    )  # Crystal receptor

    # Compute flexible residues RMSD
    flex, receptor, crystal = frmsd.load_systems(flexpath, receptorpath, crystalpath)
    r, maxr = frmsd.rmsd(flex, receptor, crystal)

    print(f"RMSD = {r:.5f}")
    print(f"MAX min(RMSD) = {maxr:.5f}")

    # Load ligand and receptor
    cmd.load(ligandpath, "ligand")  # Selection name: ligand
    cmd.load(flexpath, "flex")  # Selection name: flex
    cmd.load(crystalpath, "crystal")  # Selection name: receptor

    # Hide everything
    cmd.hide("all")

    # Show docked ligand
    cmd.show("sticks", "ligand")
    cmd.show("spheres", "ligand")
    cmd.set("sphere_scale", 0.1, "ligand")
    cmd.set("stick_radius", 0.1, "ligand")
    cmd.color("gray", "ligand")

    # Show flexible residues
    cmd.show("licorice", "flex")
    cmd.set("stick_radius", 0.3, "flex")
    cmd.color("yellow", "flex and name C*")

    # Center and zoom to ligand
    cmd.center("flex")
    cmd.zoom("flex", 5)

    # Remove solvent
    cmd.remove("solvent")

    # Get residues of receptor close to flexible residues
    stored.list = []
    cmd.iterate(
        "crystal within 0.1 of flex",  # Selection
        "stored.list.append((resn, resi, chain))",  # Action
    )

    # Remove redundancies
    flexres = set(stored.list)  # Set of flexible residues

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


cmd.extend("flexrmsd", flexrmsd)
