"""
Given a PDB file containing a list of flexible residues, get the corresponding residues
of the crystal structure.
"""

import MDAnalysis as mda
import MDAnalysis.analysis.rms as RMS

import argparse as ap
import warnings

from typing import Optional, Tuple


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

    flexres = flex.select_atoms("protein").residues

    residues = []
    for res in flexres:
        residues.append((res.resid, res.resname, res.segid))

    # Check that all flexible residues are listed
    assert len(residues) == len(flexres)

    selection = "".join(
        [f"(resid {id} and resname {name} and segid {chain}) or " for id, name, chain in residues]
    )
    selection = selection[:-4]  # Remove final " or "

    # Remove H atoms
    # TODO: Possibly need perception from atom name, when type is not present
    selection = f"({selection}) and not (type H or name H*)"

    p_atoms = protein.select_atoms(selection)
    c_atoms = crystal.select_atoms(selection)

    # for p_atom, c_atom in zip(p_atoms, c_atoms):
    #    assert p_atom.type == c_atom.type
    #    assert p_atom.name == c_atom.name
    #    assert p_atom.resid == c_atom.resid
    #    assert p_atom.resname == c_atom.resname
    #    assert p_atom.segid == c_atom.segid

    # Check that the number of atoms in the two selections is equal
    assert len(p_atoms) == len(c_atoms)

    return RMS.rmsd(p_atoms.positions, c_atoms.positions)


if __name__ == "__main__":

    args = parse()

    flex, protein, crystal = load_systems(args.flex, args.protein, args.crystal)

    r = rmsd(flex, protein, crystal)

    print(f"{r:.5f}")
