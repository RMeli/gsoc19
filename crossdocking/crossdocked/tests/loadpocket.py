from pymol import cmd, stored

import os


def loadpocket(
    pocket,
    idx,
    dataset,
    cdpath="~/Desktop/CrossDocked",
):
    # Get flexible distance from folder name
    flexdist = float(dataset.split("-")[-1][1:])

    # Clear everything
    cmd.reinitialize("everything")

    # Convert string of indices to numbers
    idx = int(idx)

    # Print informations
    print(f"Loading system {idx} in pocket {pocket}")
    print(f"flexdist = {flexdist}")

    dpath = os.path.join(dataset, "docking", pocket)

    files = [f for f in os.listdir(dpath) if f[:7] == "flexlig"]

    print(files)

    # Select ligand and receptor files
    ligname = files[idx] # Load only system IDX
    flexname = ligname.replace("flexlig", "flexrec").replace(".sdf", ".pdb")
    cligname = ligname.split("-")[1] + ".sdf"
    crecname = flexname.split("-")[2].split(".")[0] + ".pdb"

    # Build paths
    ligandpath = os.path.join(dpath, ligname)  # Docked ligand
    flexpath = os.path.join(dpath, flexname)  # Flexible residues
    cligpath = os.path.join(cdpath, pocket, cligname)
    crecpath = os.path.join(cdpath, pocket, crecname)

    # Load everything
    cmd.load(ligandpath, "ligand")  # Selection name: ligand
    cmd.load(flexpath, "flex")  # Selection name: flex
    cmd.load(crecpath, "receptor")
    cmd.load(cligpath, "cligand")

    # Hide everything
    cmd.hide("all")

    # Show receptor
    cmd.show("cartoon", "receptor")
    cmd.color("grey", "receptor")

    # Show docked ligand
    for name, color in [("ligand", "marine"), ("cligand", "pink")]:
        cmd.show("sticks", name)
        cmd.show("spheres", name)
        cmd.set("sphere_scale", 0.2, name)
        cmd.color(color, f"{name} and name C*")

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
            cmd.show("licorice", sel)
            cmd.set("stick_radius", 0.15, sel)
            cmd.color("deeppurple", sel)
            cmd.remove(f"hydro in ({sel})")

    # Show metal atoms
    cmd.show("spheres", "metals")
    cmd.set("sphere_scale", 0.5, "metals")
    
cmd.extend("loadpocket", loadpocket)
