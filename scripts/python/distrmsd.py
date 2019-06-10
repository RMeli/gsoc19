"""
Plot RMSD distributions for top 5 ligand poses across the dataset.
"""

import argparse as ap
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

from typing import Optional

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

    parser = ap.ArgumentParser(description="Compute RMSD distributions for docking poses.")

    parser.add_argument("input", type=str, help="Input file (.csv)")
    parser.add_argument("-mr", "--maxrank", type=int, default=5, help="Input file (.csv)")
    parser.add_argument("-b", "--bins", type=int, default=None, help="Input file (.csv)")
    parser.add_argument("-opath", "--outputpath", type=str, default="", help="Output files path")

    return parser.parse_args(args)


def get_lig_rmsd_for_rank(rank : int, df : pd.DataFrame) -> np.ndarray:

    return df["rmsd_lig"].loc[df["rank"] == rank].values


if __name__ == "__main__":

    import os

    args = parse()

    df = pd.read_csv(args.input)

    plt.figure()
    for rank in range(1,args.maxrank + 1):

        rmsd = get_lig_rmsd_for_rank(rank, df)

        sns.distplot(rmsd, bins=args.bins, hist=True, kde=True, label = f"rank {rank}")


    rmsd = get_lig_rmsd_for_rank(1, df)
    N = len(rmsd)
    n = len(rmsd[rmsd < 1])
    print(f"Number of sub-1Å top poses: {n} ({n / N * 100:.2f}%)")

    plt.title("RMSD Distribution")
    plt.xlabel("RMSD (Å)")
    plt.legend()
    plt.xlim([0,None])
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.pdf"))