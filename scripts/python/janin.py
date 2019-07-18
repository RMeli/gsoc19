"""
Janin plot for flexible residues.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Janin, Janin_ref

import numpy as np
from matplotlib import pyplot as plt

import argparse as ap
import os
import re
import warnings

from typing import Optional, Dict

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

def load_as_traj(system: str, rootdir: str = ""):
    """
    Load a single pose for flexible residues and all the poses in a single trajectory
    for the receptor.

    Args:
        system (str): Name of the system
        rootdir (str): Root directory for the system
    """

    # List all protein paths
    protnames = [
        os.path.join(rootdir, prot) 
        for prot in os.listdir(os.path.join(rootdir))
        if re.match(system + "_protein-?[0-9]{0,2}\.pdb", prot)
    ]

    # Load all proteins in trajectory
    prot = mda.Universe(protnames[0], protnames)

    # List all flexible residues paths
    flexnames = [
        os.path.join(rootdir, flex) 
        for flex in os.listdir(os.path.join(rootdir))
        if re.match(system + "_flex-?[0-9]{0,2}\.pdb", flex)
    ]

    # Load a single flexible residue
    flex = mda.Universe(flexnames[0])

    return flex, prot

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

def select_flexres(flex : mda.Universe, prot: mda.Universe) -> mda.AtomGroup:
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
        fres.append((res.resid, res.icode, res.segid))

    sel = "".join(
        [
            f"(resnum {id} and segid {chain}) or " 
            for id, icode, chain in fres
            if icode == ""
        ]
    ) + "".join([
            f"(resnum {id} and icode {icode} and segid {chain}) or " 
            for id, icode, chain in fres
            if icode != ""
        ]
    )

    # Sanitize selection and remove residues without Janin dihedrals
    # Ignoring them explicitly removes a warning
    sel=sel[:-4] + "and not (resname ALA or resname CYS or resname GLY or resname PRO or resname SER or resname THR or resname VAL)"
    
    return prot.select_atoms(sel)


def get_flexres_selection(
    system: str, rootdir: str = "", print_warnings : bool =False
) -> mda.AtomGroup:
    """
    For a given system, extract the flexible residues from the protein (including 
    backbone atoms) and return the corresponding selection/
    """

    if not print_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            flex, prot = load_as_traj(system, rootdir)

            flexsel = select_flexres(flex, prot)

    else:
        flex, prot = load_as_traj(system, rootdir)

        flexsel = select_flexres(flex, prot)

    return flexsel

if __name__ == "__main__":

    args = parse()

    systems = list_systems(args.systems, args.root)

    sels = [get_flexres_selection(system, path) for system, path in systems.items()]

    fig, ax = plt.subplots(figsize=plt.figaspect(1))

    X, Y = np.meshgrid(np.arange(0, 360, 6), np.arange(0, 360, 6))
    levels = [1, 6, 600]
    colors = ['#A1D4FF', '#35A1FF']
    ax.contourf(X, Y, np.load(Janin_ref), levels=levels, colors=colors)
    
    for sel in sels:
        J = Janin(sel).run()
        J.plot(ax=ax, color="k", ref=False)
    
    plt.show()