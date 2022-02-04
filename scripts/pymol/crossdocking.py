from pymol import cmd, stored

import os


def crossdocking(
    ligand,
    protein,
    pocket,
    idx="1",
    flexdist="3.5",
    path="carlos_cd",
):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idx = int(idx)

    # Convert flexdist to float
    flexdist = float(flexdist)

    # Print informations
    # print(f"Loading {dataset}/{system} (rank {idx})")
    # print(f"flexdist = {flexdist}")

    # Build paths
    ppath = os.path.join(path, pocket, "PDB_Structures")

    ligname = f"{protein}_PRO_{ligand}_LIG_aligned_v2_default_ensemble_none_flexdist3.5_p{idx}.sdf.gz"
    ligandpath = os.path.join(ppath, ligname)  # Docked ligand
    flexname = f"{protein}_PRO_{ligand}_LIG_aligned_v2_default_ensemble_none_flexdist3.5_flex_p{idx}.pdb.gz"
    flexpath = os.path.join(ppath, flexname)  # Flexible residues
    receptorname = f"{protein}_PRO_{ligand}_LIG_aligned_v2_default_ensemble_none_flexdist3.5_full_p{idx}.pdb.gz"
    receptorpath = os.path.join(ppath, receptorname)
    crystalname = f"{protein}_PRO.pdb"
    crystalpath = os.path.join(ppath, crystalname)  # Crystal receptor

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


cmd.extend("crossdocking", crossdocking)
