import MDAnalysis as mda

import argparse as ap
import numpy as np
import os
import warnings

from typing import Optional

def load_pdbqt(fname : str) -> mda.Universe:
    """
    Load .pdbqt file.

    Args:
        fname (str): Input file name (.pdbqt)

    Returns:
        Returns a `mda.Universe` containing the coordinates from `fname`.

    .. note:
        Warning for ``mda.Universe`` is ignored.
    """

    if not os.path.isfile(fname):
        raise IOError(f"{fname} does not exsist.")

    _, ext = os.path.splitext(fname)

    if ext != ".pdbqt":
        raise NameError(f"{fname} does not appear to be a .pdbqt file.")
    
    # Fixme: Redirect warning instead of suppressing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        u = mda.Universe(fname)

    return u


def ligand_center(u: mda.Universe) -> np.ndarray:
    """
    Compute center of geometry of a ligand.

    Args:
        u (mda.Universe): MDAnalysis universe

    Returns:
        A `np.ndarray` containing the coordinates of the center of geometry

    .. note:
        Assumes that ``u`` only contains the ligand (all atoms are selected).
    """
    
    ligand = u.select_atoms("all")

    return ligand.center_of_geometry()


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

    parser.add_argument("input", type=str, help="Input file (.mol2)")
    parser.add_argument("-o", "--output", type=str, help="Output file")

    return parser.parse_args(args)


if __name__ == "__main__":

    args = parse()

    u = load_pdbqt(args.input)

    center = ligand_center(u)

    center_str = f"{center[0]:.4f} {center[1]:.4f} {center[2]:.4f}"

    if args.output is None:
        print(center_str)
    else:
        with open(args.output, 'w') as fout:
            fout.write(f"{center[0]:.4f} {center[1]:.4f} {center[2]:.4f}")

