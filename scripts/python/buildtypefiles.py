"""
Build GNINA .types file (containing annotations and receptor/ligand paths) from 
.gninatypes files and a pre-computed cross-validation split.
"""

import argparse as ap
import os
import pandas as pd
import re

from typing import Optional

datasets = ["test"]
# datasets = ["refined", "other"]


def parse(args: Optional[str] = None) -> ap.Namespace:

    parser = ap.ArgumentParser()

    parser.add_argument("datapath", type=str, help="Path to database root.")
    parser.add_argument("outpath", type=str, help="Path to gninatypes root.")
    parser.add_argument("--folds", default=None, type=str, help="Folds file.")
    parser.add_argument("--min", default=2, type=float)
    parser.add_argument("--max", default=4, type=float)

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse()

    # Load folds file
    if args.folds is not None:
        df_folds = pd.read_csv(args.folds)

        n_folds = df_folds["fold"].max() + 1

        print(f"Cross-validation with {n_folds} folds.")

    # List all folders containing gninatypes files
    dirs = [
        d
        for d in os.listdir(args.outpath)
        if re.match("^....$", d) and os.path.isdir(os.path.join(args.outpath, d))
    ]

    for system in dirs:

        # fold = df_folds[df_folds["pdb"] == system]["fold"].values[0]
        # print(fold)
        # Choose at random for testing purposes:
        import numpy as np

        fold = np.random.randint(0, 3)

        # Automatically check if system is in refined or other set
        for dataset in datasets:
            scorepath = os.path.join(
                args.datapath, dataset, system, f"{system}_score.csv"
            )

            if os.path.isfile(scorepath):
                break

        df_score = pd.read_csv(scorepath)

        # Open file for correct fold
        # TODO: Open files only once
        with open(os.path.join(args.outpath, f"alltrain{fold}.types"), "a") as fout:

            # Iterate over different docking poses for a given system
            for idx, row in df_score.iterrows():

                ligname = row["name"]
                recname = ligname.replace("ligand", "protein")

                ligpath = f"{system}/{ligname}.gninatypes"
                recpath = f"{system}/{recname}.gninatypes"

                rmsd = float(row["rmsd_lig"])
                if rmsd <= args.min:
                    annotation = 0
                elif rmsd >= args.max:
                    annotation = 1
                else:  # Discard (min, max) interval
                    continue  # Skip

                score = float(row["score"])

                record = f"{annotation} {recpath} {ligpath} # {rmsd:.4f} {score:.4f}\n"

                fout.write(record)
