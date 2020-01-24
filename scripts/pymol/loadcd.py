from pymol import cmd, stored

import os


def loadcd(
    recid,
    chain,
    ligid,
    ligname,
    pocket,
    flexdist="3",
    crossdock="~/Documents/datasets/crossdocking/CrossDocked",
    cdpath="",
):

    # Clear everything
    cmd.reinitialize("everything")

    # Convert flexdist to float
    flexdist = float(flexdist)

    # Print informations
    print(f"Loading {recid}_{chain}:{ligid}_{ligname}")
    print(f"flexdist = {flexdist}")

    # Build paths
    name = f"{recid}_{chain}_rec_{ligid}_{ligname}"
    path = os.path.join(cdpath, pocket)
    ligname = f"{name}_lig.pdb"
    ligandpath = os.path.join(path, ligname)  # Docked ligand
    flexname = f"{name}_flex.pdb"
    flexpath = os.path.join(path, flexname)  # Flexible residues
    #receptorname = f"{recid}_{chain}__rec.pdb"
    #receptorpath = os.path.join(crossdock, receptorname)
    crystalname = f"{recid}_{chain}_rec.pdb"
    crystalpath = os.path.join(
        crossdock, pocket, crystalname
    )  # Crystal receptor

    # Load ligand and receptor
    cmd.load(ligandpath, "ligand")  # Selection name: ligand
    cmd.load(flexpath, "flex")  # Selection name: flex
    #cmd.load(receptorpath, "receptor")  # Selection name: receptor
    cmd.load(crystalpath, "crystal")  # Selection name: receptor

    # Hide everything
    cmd.hide("all")

    # Show receptor
    #cmd.show("cartoon", "receptor")
    #cmd.color("grey", "receptor")
    cmd.show("cartoon", "crystal")
    cmd.color("grey", "crystal")

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
    #cmd.remove("solvent")

    # Get residues of crystal close to flexible residues
    stored.list = []
    cmd.iterate(
        "crystal within 0.1 of flex",  # Selection
        "stored.list.append((resn, resi, chain))",  # Action
    )

    # Remove redundancies
    flexres = set(stored.list)  # Set of flexible residues

    # Outline flexible residues of the receptor
    #for _, resi, chain in flexres:
    #    if chain != "":  # ???
    #        sel = f"receptor and (resi {resi} in chain {chain})"
    #        cmd.show("sphere", sel)
    #        cmd.set("sphere_scale", 0.2, sel)
    #        cmd.color("yellow", sel)
    #        cmd.remove(f"hydro in ({sel})")

    # Outline flexible residues of the crystal
    for _, resi, chain in flexres:
        if chain != "":  # ???
            sel = f"crystal and (resi {resi} in chain {chain})"
            cmd.show("licorice", sel)
            cmd.set("stick_radius", 0.15, sel)
            cmd.color("grey", sel)
            cmd.remove(f"hydro in ({sel})")

    # Show metal atoms
    cmd.show("spheres", "metals")
    cmd.set("sphere_scale", 0.5, "metals")


cmd.extend("loadcd", loadcd)
