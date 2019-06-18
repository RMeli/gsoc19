import pandas as pd
import argparse as ap
import os

from typing import Optional


def parse(args: Optional[str] = None) -> ap.Namespace:

    parser = ap.ArgumentParser()

    parser.add_argument('score', type=str)
    parser.add_argument('path', type=str)
    parser.add_argument('--output', default="gnina.types")
    parser.add_argument('--min', default=2, type=float)
    parser.add_argument('--max', default=4, type=float)

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse()

    # Load score file
    csv_name = args.score
    print(f"Loading {csv_name}...", end="")
    df = pd.read_csv(csv_name)
    print("done")

    with open(args.output, "w") as fout:

        for idx, row in df.iterrows():
            if idx == 0: # Skip header
                continue

            fold = 0 # One cross validation fold

            ligname = row["name"]
            recname = ligname.replace("ligand", "protein")

            ligpath = os.path.join(args.path, f"{ligname}_{fold}.gninatypes")
            recpath = os.path.join(args.path, f"{recname}_{fold}.gninatypes")

            rmsd = float(row["rmsd_lig"])
            if rmsd <= args.min:
                annotation = 0
            elif rmsd >= args.max:
                annotation = 1
            else: # Discard (min, max) interval
                continue # Skip

            score = float(row["score"])

            record = f"{annotation} {recpath} {ligpath} # {rmsd:.4f} {score:.4f}\n"

            fout.write(record)
