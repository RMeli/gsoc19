from pymol import cmd, stored
from pymol.cgo import *

import pandas as pd

import os, sys
import re

sys.path.insert(1, "/Users/rmeli/Documents/git/rmeli/gsoc19/scripts/python/")

import inbox

# Show interpreter
# print(sys.executable)


def print_resrmsd(df, rank):

    df_rank = df[df["rank"] == rank]

    for _, row in df_rank.iterrows():
        res = row["res"]
        rmsd = row["rmsd"]

        resname, resinfo = res.split("-")

        chain = re.search("^[A-Z]", resinfo)
        resnum = re.search("[0-9]+", resinfo)
        icode = re.search("[A-Z]$", resinfo)

        print(
            f"{resname} {chain.group(0) if chain is not None else ''} {resnum.group(0):5}{icode.group(0) if icode is not None else '':1} : {rmsd:.5f} A"
        )


def draw_box(ligandpath, L):

    L2 = float(L) / 2.0  # Half box size

    lig = inbox.load(ligandpath)

    # Box center (based on ligand)
    c = inbox.center(lig)

    # Box corners
    minX = c[0] - L2
    minY = c[1] - L2
    minZ = c[2] - L2
    maxX = c[0] + L2
    maxY = c[1] + L2
    maxZ = c[2] + L2

    # Box line width
    linewidth = 2.0

    # Box color
    r, g, b = 1.0, 1.0, 1.0

    boundingBox = [
        LINEWIDTH,
        float(linewidth),
        BEGIN,
        LINES,
        COLOR,
        float(r),
        float(g),
        float(b),
        VERTEX,
        minX,
        minY,
        minZ,  # 1
        VERTEX,
        minX,
        minY,
        maxZ,  # 2
        VERTEX,
        minX,
        maxY,
        minZ,  # 3
        VERTEX,
        minX,
        maxY,
        maxZ,  # 4
        VERTEX,
        maxX,
        minY,
        minZ,  # 5
        VERTEX,
        maxX,
        minY,
        maxZ,  # 6
        VERTEX,
        maxX,
        maxY,
        minZ,  # 7
        VERTEX,
        maxX,
        maxY,
        maxZ,  # 8
        VERTEX,
        minX,
        minY,
        minZ,  # 1
        VERTEX,
        maxX,
        minY,
        minZ,  # 5
        VERTEX,
        minX,
        maxY,
        minZ,  # 3
        VERTEX,
        maxX,
        maxY,
        minZ,  # 7
        VERTEX,
        minX,
        maxY,
        maxZ,  # 4
        VERTEX,
        maxX,
        maxY,
        maxZ,  # 8
        VERTEX,
        minX,
        minY,
        maxZ,  # 2
        VERTEX,
        maxX,
        minY,
        maxZ,  # 6
        VERTEX,
        minX,
        minY,
        minZ,  # 1
        VERTEX,
        minX,
        maxY,
        minZ,  # 3
        VERTEX,
        maxX,
        minY,
        minZ,  # 5
        VERTEX,
        maxX,
        maxY,
        minZ,  # 7
        VERTEX,
        minX,
        minY,
        maxZ,  # 2
        VERTEX,
        minX,
        maxY,
        maxZ,  # 4
        VERTEX,
        maxX,
        minY,
        maxZ,  # 6
        VERTEX,
        maxX,
        maxY,
        maxZ,  # 8
        END,
    ]

    cmd.load_cgo(boundingBox, "box")


def flexrmsd(
    system,
    dataset,
    idx="1",
    box_size="23.5",
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
    df_rmsd = pd.read_csv(os.path.join(dataset, system, f"{system}_score.csv"))
    df_recrmsd = pd.read_csv(os.path.join(dataset, system, "resrmsd.csv"))

    print_resrmsd(df_recrmsd, idx)

    print(df_rmsd[df_rmsd["rank"] == idx])

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

    # Draw box
    draw_box(ligandpath, box_size)


cmd.extend("flexrmsd", flexrmsd)
