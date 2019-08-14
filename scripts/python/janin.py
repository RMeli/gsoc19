"""
Janin plot for flexible residues.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Janin, Janin_ref

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm

import argparse as ap
import os
import re
import warnings

from typing import Optional, Dict

# Colormap for different ranks
cmap = cm.get_cmap("magma")
colors = {i + 1: cmap(c) for i, c in enumerate(np.linspace(0, 1, 20))}


def parse(args: Optional[str] = None) -> ap.Namespace:
    """
    Parsed command line arguments.

    Args:
        args (str, optional): String to parse

    Returns:
        An `ap.Namespace` containing the parsed options

    .. note::
        If ``args is None``, parse from ``sys.argv``
    """

    parser = ap.ArgumentParser(description="Compute center of geometry of a molecule.")

    parser.add_argument("systems", type=str, help="List of paths to systems")
    parser.add_argument("-r", "--root", type=str, default="", help="")

    return parser.parse_args(args)


def list_poses(system: str, rootdir: str = ""):

    # List all protein paths
    protnames = [
        os.path.join(rootdir, prot)
        for prot in os.listdir(os.path.join(rootdir))
        if re.match(system + "_protein-[0-9]{1,2}\.pdb", prot)
    ]

    # List all flexible residues paths
    flexnames = [name.replace("protein", "flex") for name in protnames]

    # Extract ranks
    names = [os.path.basename(name).replace(".pdb", "") for name in flexnames]
    ranks = [int(name.replace(f"{system}_flex-", "")) for name in names]

    return zip(ranks, flexnames, protnames)


def list_systems(syslist: str, rootdir: str = "") -> Dict[str, str]:
    """
    List all systems in a given list and store them in a dictionary (as key) along with
    the path to the system (as value).

    Args:
        syslist (str): Name of the files listing the systems
        rootdir (str): Root directory (to be appended to the path in `syslist`)

    Returns:
        Dictionaty containing the systems' names as keys and corresponding paths as
        valuse
    """

    systems = {}
    with open(syslist, "r") as file:
        for line in file:
            l = line.strip()
            system = os.path.basename(l)
            systems[system] = os.path.join(rootdir, l)

    return systems


def select_flexres(flex: mda.Universe, prot: mda.Universe) -> mda.AtomGroup:
    """
    Given a protein and a series of flexible residues, selectss the full flexible
    residues (including backbone atoms) from the protein structure.

    Args:
        flex (mda.Universe): flexible residues
        prot (mda.Universe): protein

    Returns:
        An `mda.AtomGroup` containing the atoms corresponding to flexible residues
        extracted from the protein (including backbone atoms)
    """

    fres = []
    for res in flex.residues:
        fres.append((res.resnum, res.resname, res.icode, res.segid))

    sel = "".join(
        [
            f"(resid {num}{icode} and resname {name} and segid {chain}) or "
            for num, name, icode, chain in fres
        ]
    )

    # Sanitize selection and remove residues without Janin dihedrals
    # Ignoring them explicitly removes a warning
    sel = (
        sel[:-4]
        + "and not (resname ALA or resname CYS or resname GLY or resname PRO or resname SER or resname THR or resname VAL)"
    )

    return prot.select_atoms(sel)


def plot_reference(ax, ref=None):
    if ref is None:
        ref = np.load(Janin_ref)

    X, Y = np.meshgrid(np.arange(0, 360, 6), np.arange(0, 360, 6))
    levels = [1, 6, 600]
    colors = ["#A1D4FF", "#35A1FF"]
    ax.contourf(X, Y, ref, levels=levels, colors=colors)


def plot(ax, systems):

    for system, path in systems.items():
        for rank, flexname, protname in list_poses(system, path):
            print(flexname, protname)

            flex = mda.Universe(flexname)
            prot = mda.Universe(protname)

            sel = select_flexres(flex, prot)

            J = Janin(sel).run()
            J.plot(ax=ax, color=colors[rank], ref=False)


if __name__ == "__main__":

    args = parse()

    systems = list_systems(args.systems, args.root)

    # sels = [get_flexres_selection(system, path) for system, path in systems.items()]

    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    plot_reference(ax)

    plot(ax, systems)

    cm = matplotlib.colors.LinearSegmentedColormap.from_list(
        "cm", list(colors.values())
    )
    Z = [[0, 0], [0, 0]]
    levels = range(1, 22)
    im = ax.contourf(Z, levels, cmap=cm)
    fig.colorbar(im)

    plt.show()
