"""
Build GNINA .types file (containing annotations and receptor/ligand paths) from 
.gninatypes files for the different systems.
"""

from sklearn.model_selection import KFold

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
    parser.add_argument("--min", default=2, type=float)
    parser.add_argument("--max", default=4, type=float)

    args = parser.parse_args()

    return args


def write_record(df_score: pd.DataFrame, min: float, max: float, outfile):

    # Iterate over different docking poses for a given system
    for _, row in df_score.iterrows():

        ligname = f"{row['system']}_ligand-{row['rank']}"
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

            write_record(df_score, args.min, args.max, out) # Write record