"""
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
    parser.add_argument("-opath", "--outputpath", type=str, default="", help="Output files path")

    return parser.parse_args(args)


def get_lig_rmsd_for_rank(rank : int, df : pd.DataFrame) -> np.ndarray:

    return df["rmsd_lig"].loc[df["rank"] == rank].values


if __name__ == "__main__":

    import os

    args = parse()

    df = pd.read_csv(args.input)

    plt.figure()
    for rank in range(1,3):
        name = os.path.join(args.outputpath, f"distrmsd-{rank}.pdf")

        rmsd = get_lig_rmsd_for_rank(rank, df)

        sns.distplot(rmsd, hist=True, kde=True, label = f"rank {rank}")

    plt.title("RMSD Distribution")
    plt.xlabel("RMSD (Å)")
    plt.legend()
    plt.xlim([0,None])
    plt.show()