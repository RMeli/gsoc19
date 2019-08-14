"""
Build GNINA .types file (containing annotations and receptor/ligand paths) from 
.gninatypes files for the different systems.
"""

import argparse as ap
import os
import pandas as pd
import re

from typing import Optional

datasets = ["refined", "other"]


def parse(args: Optional[str] = None) -> ap.Namespace:

    parser = ap.ArgumentParser()

    parser.add_argument("datapath", type=str, help="Path to database root.")
    parser.add_argument("outpath", type=str, help="Path to gninatypes root.")
    parser.add_argument("--lmin", default=2, type=float)
    parser.add_argument("--lmax", default=4, type=float)
    parser.add_argument("--fmin", default=1, type=float)
    parser.add_argument("--fmax", default=1.5, type=float)

    args = parser.parse_args()

    return args


def write_record(
    system: str,
    df_score: pd.DataFrame,
    lmin: float,
    lmax: float,
    fmin: float,
    fmax: float,
    outfile: str,
):

    # Iterate over different docking poses for a given system
    for _, row in df_score.iterrows():

        rank = int(round(row["rank"]))

        ligname = f"{system}_ligand-{rank}"
        recname = ligname.replace("ligand", "protein")

        ligpath = f"{system}/{ligname}.gninatypes"
        recpath = f"{system}/{recname}.gninatypes"

        rmsd_lig = float(row["rmsd_lig"])
        rmsd_flex = float(row["rmsd_fmax"])
        if rmsd_lig <= lmin and rmsd_flex <= fmin:
            annotation = 0
        elif rmsd_lig >= lmax and rmsd_flex >= fmax:
            annotation = 1
        else:  # Discard (lmin, lmax) and (fmin, fmax) intervals
            continue  # Skip

        score = float(row["score"])

        record = f"{annotation} {recpath} {ligpath} # {rmsd_lig:.4f}  {rmsd_flex:.4f} {score:.4f}\n"

        outfile.write(record)


if __name__ == "__main__":

    args = parse()

    # List all folders containing gninatypes files
    dirs = [
        d
        for d in os.listdir(args.outpath)
        if re.match("^....$", d) and os.path.isdir(os.path.join(args.outpath, d))
    ]

    # Loop over folds
    for system in dirs:

        typesfile = f"{args.outpath}/{system}/{system}.types"

        with open(typesfile, "w") as out:

            # Automatically check if system is in refined or other set
            for dataset in datasets:
                scorepath = os.path.join(
                    args.datapath, dataset, system, f"{system}_score.csv"
                )

                if os.path.isfile(scorepath):
                    break

            # Get system scores
            df_score = pd.read_csv(scorepath)

            write_record(
                system, df_score, args.lmin, args.lmax, args.fmin, args.fmax, out
            )  # Write record
