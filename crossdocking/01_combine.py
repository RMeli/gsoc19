import numpy as np
import os

from collections import namedtuple, defaultdict
from typing import List

Receptor = namedtuple("Receptor", ["pdbid", "chain"])
Ligand = namedtuple("Ligand", ["pdbid", "name"])

def receptorsdict(rfname: str):

    receptors = defaultdict(list)

    with open(rfname, "r") as fin:
        for line in fin:
            pocket, pdbid, chain = line.strip().split(":")

            receptors[pocket].append(Receptor(pdbid=pdbid, chain=chain))

    return receptors


def ligandsdict(lfname: str):

    ligands = defaultdict(list)

    with open(lfname, "r") as fin:
        for line in fin:
            pocket, pdbid, name, _ = line.strip().split(":")

            ligands[pocket].append(Ligand(pdbid=pdbid, name=name))

    return ligands


def ligandpath(lig, pocket, root):

    ligprefix = f"{lig.pdbid}_{lig.name}"

    ligpath = ""
    try:
        ligname = ligprefix + "_uff2.sdf"
        ligpath = os.path.join(root, pocket, ligname)

        open(ligpath, "r")
    except FileNotFoundError:
        ligname = ligprefix + "_uff.sdf"
        ligpath = os.path.join(root, pocket, ligname)

    return ligpath



def crossdocking(ligdict, recdict, outfile: str = "crossdocking.dat", root: str = ""):

    with open(outfile, "w") as fout:
        for pocket, reclist in recdict.items():
            for rec in reclist:

                ligands = ligdict[pocket]

                for idx, lig in enumerate(ligdict[pocket]):

                    recname = f"{rec.pdbid}_{rec.chain}_rec.pdb"
                    recpath = os.path.join(root, pocket, recname)

                    ligpath = ligandpath(lig, pocket, root)

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

        return parser.parse_args(args)

    args = parse()

    ligands = ligandsdict(args.lfile)
    receptors = receptorsdict(args.rfile)

    crossdocking(ligands, receptors, args.output, args.root)
