"""
Build GNINA .types file (containing annotations and receptor/ligand paths) from 
.gninatypes files and a pre-computed cross-validation split.
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
    parser.add_argument("--folds", default=None, type=str, help="Folds file.")
    parser.add_argument("--min", default=2, type=float)
    parser.add_argument("--max", default=4, type=float)

    args = parser.parse_args()

    return args


def write_record(df_score, outfile):

    # Iterate over different docking poses for a given system
    for _, row in df_score.iterrows():

        ligname = f"{row['system']}_ligand-{row['rank']}.pdb"
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

    # Load folds file
    if args.folds is not None:
        df_folds = pd.read_csv(args.folds)

        n_folds = df_folds["fold"].max() + 1

        print(f"Cross-validation with {n_folds} folds.")

    else:
        raise Exception # TODO: Support random splitting

    # List all folders containing gninatypes files
    dirs = [
        d
        for d in os.listdir(args.outpath)
        if re.match("^....$", d) and os.path.isdir(os.path.join(args.outpath, d))
    ]

    kf = KFold(n_splits=n_folds)

    # Loop over folds
    for train, test in kf.split(range(n_folds), range(n_folds)):
        idx = test[0]

        print(f"   Building train {train} and test {test}...")

        # Open type files for current fold
        trainfile = f"{args.outpath}/alltrain{idx}.types"
        testfile = f"{args.outpath}/alltest{idx}.types"
        with open(trainfile, "w") as trainout, open(testfile, "w") as testout:

            # Loop over systems
            for system in dirs:

                # Automatically check if system is in refined or other set
                for dataset in datasets:
                    scorepath = os.path.join(
                        args.datapath, dataset, system, f"{system}_score.csv"
                    )

                    if os.path.isfile(scorepath):
                        break

                # Get system scores
                df_score = pd.read_csv(scorepath)

                # Get actual fold of the system (pre-computed)
                try:
                    fold = df_folds[df_folds["pdb"] == system]["fold"].values[0]
                except IndexError:
                    print(f"WARNING: Skypping system {system} with unknown fold...")
                    continue

                if fold in train: # System should be in train
                    write_record(df_score, trainout) # Write on train file
                elif fold in test: # System should be in test
                    write_record(df_score, testout) # Write on test file
                else:
                    raise ValueError(f"Fold {fold} is invalid.")