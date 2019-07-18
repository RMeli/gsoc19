import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Janin

from matplotlib import pyplot as plt

import argparse as ap
import os
import re
import warnings

from typing import Optional

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

    parser.add_argument("system", type=str, help="")
    parser.add_argument("-r", "--root", type=str, default="", help="")

    return parser.parse_args(args)

def load_as_traj(system: str, rootdir: str = ""):

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

def select_flexres(flex : mda.Universe, prot: mda.Universe) -> mda.AtomGroup:

    fres = []
    for res in flex.residues:
        fres.append((res.resid, res.icode, res.segid))

    print(fres)

    sel = "".join(
        [
            f"(resnum {id} and segid {chain}) or " 
            for id, icode, chain in fres
            if icode == ""
        ]
    ).join(
        [
            f"(resnum {id} and icode {icode} and segid {chain}) or " 
            for id, icode, chain in fres
            if icode != ""
        ]
    )

    sel=sel[:-4]

    print(sel)

    return prot.select_atoms(sel)


def get_flexres_selection(
    system: str, rootdir: str = "", print_warnings : bool =False
) -> mda.Universe:

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

    sel = get_flexres_selection(args.system, args.root)
    
    J = Janin(sel).run()

    fig, ax = plt.subplots(figsize=plt.figaspect(1))
    J.plot(ax=ax, color='k', ref=True)
    plt.show()

