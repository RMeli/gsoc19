import MDAnalysis as mda # MDA 2.0 or higher

from rdkit import Chem

import argparse as ap
import numpy as np
import os
import warnings
import gzip

from typing import Optional, Tuple

def loadsdf(fname: str, gzipped=False) -> mda.Universe:
    if not gzipped:
            rdmol = next(Chem.SDMolSupplier(fname, removeHs=False))
    else:
        with gzip.open(fname, "r") as fgz:
            rdmol = next(Chem.ForwardSDMolSupplier(fgz, removeHs=False))

    # Convert RDKit molecule to MDA universe
    return mda.Universe(rdmol)

def loadpdb(fname: str) -> mda.Universe:

    if not os.path.isfile(fname):
        raise IOError(f"{fname} does not exsist.")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        u = mda.Universe(fname)

    return u


def center(u: mda.Universe) -> np.ndarray:

    return u.atoms.center_of_geometry()


def min_xyz(u: mda.Universe) -> np.ndarray:

    coords = u.trajectory[0].positions

    return np.min(coords, axis=0)


def max_xyz(u: mda.Universe) -> np.ndarray:

    coords = u.trajectory[0].positions

    return np.max(coords, axis=0)


def in_box(c, mmin, mmax, L) -> bool:

    L2 = L / 2.0  # Half box size

    dplus = mmax - c
    dminus = c - mmin

    return True if np.alltrue(dplus < L2) and np.alltrue(dminus < L2) else False

def inbox(ligname: str, flexname: str, box_size: float) -> Tuple[bool, bool]:

    lignoext, ligext = os.path.splitext(ligname)
    if ligext == ".gz":
        gzipped=True
        _, ligext = os.path.splitext(lignoext)
    else:
        gzipped=False

    if ligext == ".pdb":
        ligand = loadpdb(ligname)
    elif ligext == ".sdf":
        ligand = loadsdf(ligname, gzipped)
    else:
        raise Exception(f"Extension {ligext} not supported.")

    flex = loadpdb(flexname)

    # Ligand center (box center)
    c = center(ligand)

    # Ligand min and max positions
    lmin = min_xyz(ligand)
    lmax = max_xyz(ligand)

    # Flexible residues min and max positions
    fmin = min_xyz(flex)
    fmax = max_xyz(flex)

    ligin: bool = in_box(c, lmin, lmax, box_size)
    flexin: bool = in_box(c, fmin, fmax, box_size)

    return ligin, flexin

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

    parser.add_argument("ligand", type=str, help="Ligand file (.mol2)")
    parser.add_argument("flex", type=str, help="Flexible residues file (PDB)")
    parser.add_argument("-L", "--box_size", type=float, default=23.5, help="Box size")

    return parser.parse_args(args)


if __name__ == "__main__":

    args = parse()

    ligin, flexin = inbox(args.ligand, args.flex, args.box_size)

    print(ligin, flexin)