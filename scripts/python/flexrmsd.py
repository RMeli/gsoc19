"""
Given a PDB file containing a list of flexible residues, get the corresponding residues
of the crystal structure.
"""

import MDAnalysis as mda
import MDAnalysis.analysis.rms as RMS

import argparse as ap
import warnings

from typing import Optional, Tuple

import sys


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

    parser.add_argument("flex", type=str, help="Flexible residues (PDB)")
    parser.add_argument("protein", type=str, help="Protein structure (PDB)")
    parser.add_argument("crystal", type=str, help="Crystal structure (PDB)")
    parser.add_argument(
        "-o", "--output", type=str, default="cflex.pdb", help="Output residues"
    )

    return parser.parse_args(args)


def load_systems(
    flexname: str, proteinname: str, crystalname: str,
    print_warnings=False
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


def rmsd(flex: mda.Universe, protein: mda.Universe, crystal: mda.Universe) -> float:

    def sel(resnum, resname, segid, icode) -> str:
        s = f"(resid {resnum}{icode} and resname {resname} and segid {segid})"

        return s

    flexres = flex.select_atoms("protein").residues

    max_rmsd = -1

    residues = []
    selection = ""
    for res in flexres:
        ressel = sel(res.resnum, res.resname, res.segid, res.icode) + " and not (type H or name H*)"

        # Select single residue
        p_res = protein.select_atoms(ressel)
        c_res = crystal.select_atoms(ressel)

        # Compute minimised RMSD for single residue
        res_rmsd = RMS.rmsd(p_res.positions, c_res.positions, superposition=True)
            
        # Store the maximum RMSD for a single residue
        if res_rmsd > max_rmsd:
            max_rmsd = res_rmsd

        residues.append((res.resnum, res.resname, res.segid, res.icode))

    # Check that all flexible residues are listed
    assert len(residues) == len(flexres)

    # TODO: Can be improved by using ressel
    selection = "".join(
        [ 
            sel(id, name, chain, icode) + " or "
            for id, name, chain, icode in residues
        ]
    )
    selection = selection[:-4]  # Remove final " or "

    # Remove H atoms
    # TODO: Possibly need perception from atom name, when type is not present
    selection = f"({selection}) and not (type H or name H*)"

    p_atoms = protein.select_atoms(selection)
    c_atoms = crystal.select_atoms(selection)

    # Check that the number of atoms in the two selections is equal
    assert len(p_atoms) == len(c_atoms)

    return RMS.rmsd(p_atoms.positions, c_atoms.positions), max_rmsd

if __name__ == "__main__":

    args = parse()

    flex, protein, crystal = load_systems(args.flex, args.protein, args.crystal)

    r, maxr = rmsd(flex, protein, crystal)

    print(f"{r:.5f} {maxr:.5}")
