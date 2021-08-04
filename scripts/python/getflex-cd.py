"""
Split flexible residues into different files to compute per-residue RMSD>
"""

import MDAnalysis as mda

import argparse as ap
import warnings

from typing import Optional, Tuple

import sys, os


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

    parser = ap.ArgumentParser(description="Get flexible residues.")

    parser.add_argument("flex", type=str, help="Flexible residues (PDB)")
    parser.add_argument("cflex", type=str, help="Crystal structure (PDB)")
    parser.add_argument("-d", "--dir", type=str, default="", help="Root directory")
    parser.add_argument("-r", "--rank", type=str, default=None, help="Rank")

    return parser.parse_args(args)


def load_systems(
    flexname: str, cflexname: str, print_warnings=False
) -> Tuple[mda.Universe, mda.Universe, mda.Universe]:

    if not print_warnings:
        with warnings.catch_warnings():
            # TODO: log warnings instead of ignoring
            warnings.simplefilter("ignore")

            flex = mda.Universe(flexname)
            cflex = mda.Universe(cflexname)
    else:
        flex = mda.Universe(flexname)
        cflex = mda.Universe(cflexname)

    # Add cell dimensions to avoid warnings
    flex.dimensions=[1,1,1,90,90,90]
    cflex.dimensions=[1,1,1,90,90,90]

    return flex, cflex


def save_systems(
    flex: mda.Universe, cflex: mda.Universe, dir: str, rank: str,
):
    def sel(resnum, resname, segid, icode) -> str:
        s = f"(resid {resnum}{icode} and resname {resname} and segid {segid})"

        return s

    flexres = flex.select_atoms("protein").residues

    residues = []
    for res in flexres:
        ressel = (
            sel(res.resnum, res.resname, res.segid, res.icode)
            + " and not (type H or name H*)"
        )

        # Select single residue
        p_res = flex.select_atoms(ressel)
        c_res = cflex.select_atoms(ressel)

        assert p_res.n_atoms == c_res.n_atoms

        if rank != "":
            srank=f"-p{rank}"
        else:
            srank = ""

        pfname = os.path.join(
            dir, f"pflex{srank}-{res.resname}-{res.segid}{res.resnum}{res.icode}.pdb"
        )
        cfname = os.path.join(
            dir, f"cflex{srank}-{res.resname}-{res.segid}{res.resnum}{res.icode}.pdb"
        )

        # Write out PDB files
        p_res.write(pfname)
        c_res.write(cfname)

if __name__ == "__main__":

    args = parse()

    flex, cflex = load_systems(args.flex, args.cflex)

    save_systems(flex, cflex, args.dir, args.rank)
