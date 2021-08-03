"""
Given a PDB file containing a list of flexible residues, get the corresponding residues
of the crystal structure.
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
    parser.add_argument("protein", type=str, help="Protein structure (PDB)")
    parser.add_argument("crystal", type=str, help="Crystal structure (PDB)")
    parser.add_argument("-d", "--dir", type=str, default="", help="Root directory")

    return parser.parse_args(args)


def load_systems(
    flexname: str, proteinname: str, crystalname: str, print_warnings=False
) -> Tuple[mda.Universe, mda.Universe, mda.Universe]:

    if not print_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            flex = mda.Universe(flexname)
            protein = mda.Universe(proteinname)
            crystal = mda.Universe(crystalname)

    else:
        flex = mda.Universe(flexname)
        protein = mda.Universe(proteinname)
        crystal = mda.Universe(crystalname)

    return flex, protein, crystal


def save_systems(
    flex: mda.Universe, protein: mda.Universe, crystal: mda.Universe, dir: str
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
        p_res = protein.select_atoms(ressel)
        c_res = crystal.select_atoms(ressel)

        assert p_res.n_atoms == c_res.n_atoms

        pfname = os.path.join(
            dir, f"pflex-{res.resname}-{res.segid}{res.resnum}{res.icode}.pdb"
        )
        cfname = os.path.join(
            dir, f"cflex-{res.resname}-{res.segid}{res.resnum}{res.icode}.pdb"
        )

        # Write out PDB files
        p_res.write(pfname)
        c_res.write(cfname)

        residues.append((res.resnum, res.resname, res.segid, res.icode))

    # Check that all flexible residues are listed
    assert len(residues) == len(flexres)

    # TODO: Can be improved by using ressel
    selection = "".join(
        [sel(id, name, chain, icode) + " or " for id, name, chain, icode in residues]
    )
    selection = selection[:-4]  # Remove final " or "

    # Remove H atoms
    # TODO: Possibly need perception for atom name, when type is not present
    selection = f"({selection}) and not (type H or name H*)"

    p_atoms = protein.select_atoms(selection)
    c_atoms = crystal.select_atoms(selection)

    # Check that the number of atoms in the two selections is equal
    assert len(p_atoms) == len(c_atoms)

    pfname = os.path.join(dir, "pflex.pdb")
    cfname = os.path.join(dir, "cflex.pdb")

    p_atoms.write(pfname)
    c_atoms.write(cfname)


if __name__ == "__main__":

    args = parse()

    flex, protein, crystal = load_systems(args.flex, args.protein, args.crystal)

    save_systems(flex, protein, crystal, args.dir)
