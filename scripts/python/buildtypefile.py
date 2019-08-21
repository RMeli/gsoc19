"""
Build GNINA .types file (containing annotations and receptor/ligand paths) from 
.gninatypes files for the different systems.
"""

import argparse as ap
import os
import pandas as pd
import re

from inbox import inbox

from typing import Optional

def parse(args: Optional[str] = None) -> ap.Namespace:

    parser = ap.ArgumentParser()

    parser.add_argument("datapath", type=str, help="Path to database root.")
    parser.add_argument("typespath", type=str, help="Path to gninatypes root.")
    parser.add_argument("--lmin", default=2, type=float)
    parser.add_argument("--lmax", default=4, type=float)
    parser.add_argument("--fmin", default=1, type=float)
    parser.add_argument("--fmax", default=1.5, type=float)
    parser.add_argument("-L", "--box_size", default=None, type=float)
    parser.add_argument("-o", "--out", default="all.types", type=str)
    parser.add_argument("-d", "--datasets", nargs="+", default=["refined", "other"], type=str)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("--all", action="store_true", default=False, help="Do not exclude poses in [lmin,lmax] and [fmin,fmax]")

    args = parser.parse_args()

    return args


def write_record(
    system: str,
    df_score: pd.DataFrame,
    lmin: float,
    lmax: float,
    fmin: float,
    fmax: float,
    keepall: bool,
    box_size: float,
    path: float,
    outfile: str,
):

    # Iterate over different docking poses for a given system
    for _, row in df_score.iterrows():

        rank = int(round(row["rank"]))

        ligname = f"{system}_ligand-{rank}"
        recname = ligname.replace("ligand", "protein")
        flexname = ligname.replace("ligand", "flex")

        ligpath = f"{system}/{ligname}.gninatypes"
        recpath = f"{system}/{recname}.gninatypes"

        ligpdbpath = os.path.join(path, f"{system}_ligand-{rank}.pdb")
        flexpdbpath = os.path.join(path, f"{system}_flex-{rank}.pdb")

        rmsd_lig = float(row["rmsd_lig"])
        rmsd_flex = float(row["rmsd_fmax"])
        if rmsd_lig <= lmin and rmsd_flex <= fmin:
            annotation = 1 # Positive label

            # Check that ligand and flexible residues are within the box
            if not keepall:
                ligin, flexin = inbox(ligpdbpath, flexpdbpath, box_size)
                if ligin == False or flexin == False:
                    continue

        elif rmsd_lig >= lmax and rmsd_flex >= fmax:
            annotation = 0 # Negative label

            # Check that ligand and flexible residues are within the box
            if not keepall:
                ligin, flexin = inbox(ligpdbpath, flexpdbpath, box_size)
                if ligin == False or flexin == False:
                    continue

        else:  # Discard (lmin, lmax) and (fmin, fmax) intervals
            if not keepall:
                continue  # Skip
            else:
                annotation = 0

        score = float(row["score"])

        record = f"{annotation} {recpath} {ligpath} # {rmsd_lig:.4f} {rmsd_flex:.4f} {score:.4f}\n"

        outfile.write(record)


if __name__ == "__main__":

    args = parse()

    if args.verbose:
        print("BUILDTYPEFILE")
        print(f"  datapath =", args.datapath)
        print(f"  typespath =", args.typespath)
        print(f"  datasets =", args.datasets)
        print(f"  lmin = {args.lmin:.2f}\tlmax = {args.lmax:.2f}")
        print(f"  fmin = {args.fmin:.2f}\tfmax = {args.fmax:.2f}")
        print(f"  box_size = {args.box_size:.2f}")
        if args.all:
            print("  WARNING: Including all systems in the typefile.")

    # List all folders containing gninatypes files
    dirs = [
        d
        for d in os.listdir(args.typespath)
        if re.match("^....$", d) and os.path.isdir(os.path.join(args.typespath, d))
    ]

    with open(args.out, "w") as out:

        # Loop over folds
        for system in dirs:

            # Automatically check if system is in refined or other set
            for dataset in args.datasets:
                scorepath = os.path.join(
                    args.datapath, dataset, system, f"{system}_score.csv"
                )

                if os.path.isfile(scorepath):
                    break

            # Get system scores
            df_score = pd.read_csv(scorepath)

            # Get  path to database for ligand and flex files
            path = os.path.join(args.datapath, dataset, system)

            write_record(
                system, df_score, args.lmin, args.lmax, args.fmin, args.fmax, args.all, args.box_size, path, out
            )  # Write record
