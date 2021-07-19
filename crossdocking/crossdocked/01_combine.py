"""
Combine every ligand associated to a pocket with the cognate receptor and a given
number of random receptors with the same pocket.
"""

import numpy as np
import os

from collections import namedtuple, defaultdict
from typing import List

Receptor = namedtuple("Receptor", ["pdbid", "chain"])
Ligand = namedtuple("Ligand", ["pdbid", "name"])


def receptorsdict(rfname: str):
    """
    Get all receptors from file.

    Receptors are stored in a dictionary of lists, keyed by pocket.
    """

    receptors = defaultdict(list)

    with open(rfname, "r") as fin:
        for line in fin:
            pocket, pdbid, chain = line.strip().split(":")

            receptors[pocket].append(Receptor(pdbid=pdbid, chain=chain))

    return receptors


def ligandsdict(lfname: str):
    """
    Get all ligands from file.

    Ligands are stored in a dictionary of lists, keyed by pocket.
    """

    ligands = defaultdict(list)

    with open(lfname, "r") as fin:
        for line in fin:
            pocket, pdbid, name, _ = line.strip().split(":")

            ligands[pocket].append(Ligand(pdbid=pdbid, name=name))

    return ligands


def ligandpath(lig, pocket, root):

    ligprefix = f"{lig.pdbid}_{lig.name}"

    ligfound = False
    for ext in ["uff2", "uff", "lig_prody", "lig_babel"]:
        ligname = ligprefix + f"_{ext}.sdf"
        ligpath = os.path.join(pocket, ligname)

        if os.path.isfile(os.path.join(root, ligpath)):
            ligfound = True
            break

    if not ligfound:
        raise FileNotFoundError(
            f"\n\tpocket: {pocket} | pdbid : {lig.pdbid} | name: {lig.name}"
        )

    return ligpath


def crossdocking(
    ligdict, recdict, outfile: str = "crossdocking.dat", root: str = "", nmax: int = 2
):

    with open(outfile, "w") as fout:

        # Loop over all pockets
        for pocket, liglist in ligdict.items():

            # Loop over all ligands within a pocket
            for lig in liglist:

                ligpath = ligandpath(lig, pocket, root)

                # Get all receptors
                reclist = recdict[pocket]

                # Number of receptors
                n = len(reclist)

                # Receptor indices
                idxs = np.arange(n)

                # Find cognate receptor
                idx_cognate = np.nan
                for idx, rec in enumerate(reclist):
                    if lig.pdbid == rec.pdbid:
                        # Store cognate receptor idx for deletion
                        idx_cognate = idx

                        # Write ligand and cognate receptor
                        # This corresponds to re-docking
                        recname = f"{rec.pdbid}_{rec.chain}_rec.pdb"
                        recpath = os.path.join(pocket, recname)
                        fout.write(f"{ligpath} {recpath}\n")

                        break

                # Delete cognate receptor index from indices to be sampled
                idxs = np.delete(idxs, idx_cognate)

                # Sample nmax-1 random receptors
                samples = np.random.choice(
                    idxs, size=nmax - 1 if nmax - 1 < n - 1 else n - 1, replace=False
                )

                for idx in samples:

                    rec = reclist[idx]

                    recname = f"{rec.pdbid}_{rec.chain}_rec.pdb"
                    recpath = os.path.join(pocket, recname)

                    fout.write(f"{ligpath} {recpath}\n")


if __name__ == "__main__":

    import argparse as ap
    from typing import Optional

    def parse(args: Optional[str] = None) -> ap.Namespace:

        parser = ap.ArgumentParser(
            description="Combine receptors and ligands for cross-docking"
        )

        parser.add_argument(
            "lfile",
            type=str,
            help="Ligands file (<pocket>_0:<ligand pdb>:<3 letter code for lig in pdb>:)",
        )
        parser.add_argument(
            "rfile", type=str, help="Receptors file (<pocket>_0:<receptor pdb>:<chain>)"
        )
        parser.add_argument("-r", "--root", type=str, default="", help="Data root")
        parser.add_argument(
            "-o", "--output", type=str, default="crossdocking.dat", help=""
        )

        parser.add_argument(
            "-m",
            "--max",
            type=int,
            default=2,
            help="Maximum number of ligand-receptor pairs",
        )

        return parser.parse_args(args)

    args = parse()

    ligands = ligandsdict(args.lfile)
    receptors = receptorsdict(args.rfile)

    crossdocking(ligands, receptors, args.output, args.root, args.max)
