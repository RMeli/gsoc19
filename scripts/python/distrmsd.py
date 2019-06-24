"""
Plot RMSD distributions for top ligand poses across the dataset, given a CSV file
containing the RMSDs for all systems and all ranks.
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

    parser = ap.ArgumentParser(
        description="Compute RMSD distributions for docking poses."
    )

    parser.add_argument("input", type=str, help="Input file (.csv)")
    parser.add_argument(
        "-mr", "--maxrank", type=int, default=5, help="Input file (.csv)"
    )
    parser.add_argument(
        "-b", "--bins", type=int, default=None, help="Input file (.csv)"
    )
    parser.add_argument(
        "-opath", "--outputpath", type=str, default="", help="Output files path"
    )

    return parser.parse_args(args)


def get_lig_rmsd_for_rank(rank: int, rmsd_name: str, df: pd.DataFrame) -> np.ndarray:

    return df[rmsd_name].loc[df["rank"] == rank].values


def plot(df: pd.DataFrame, rmsd_name: str, maxrank: int, bins: int, ax):

    for rank in range(1, maxrank + 1):

        rmsd = get_lig_rmsd_for_rank(rank, rmsd_name, df)

        sns.distplot(rmsd, bins=bins, hist=True, kde=True, label=f"rank {rank}")


if __name__ == "__main__":

    import os

    args = parse()

    df = pd.read_csv(args.input)

    bad_values = [np.nan, np.inf, -np.inf]
    if df.isin(bad_values).sum(axis=0).max() > 0:
        idx = df.index[df["rmsd_lig"].isin(bad_values)]  # Get row indices of bad values

        print(f"WARNING: {args.input} contains NaN or Inf:")
        print(df.iloc[idx])

        # Replace infinite values with NaNs and drop them
        df = df.replace([np.inf, -np.inf], np.nan)  # Change Inf with NaN
        df = df.dropna()  # Drop rows with NaN

    plt.figure(figsize=(15, 8))

    ax = plt.subplot(1, 3, 1)
    plot(df, "rmsd_lig", args.maxrank, args.bins, ax)
    ax.set_xlabel("RMSD (Å)")
    ax.set_xlim([0, None])
    ax.title.set_text("Ligand")
    ax.legend()

    ax = plt.subplot(1, 3, 2)
    plot(df, "rmsd_flex", args.maxrank, args.bins, ax)
    ax.set_xlabel("RMSD (Å)")
    ax.set_xlim([0, None])
    ax.title.set_text("Flexible Residues")
    ax.legend()

    ax = plt.subplot(1, 3, 3)
    plot(df, "rmsd_tot", args.maxrank, args.bins, ax)
    ax.set_xlabel("RMSD (Å)")
    ax.set_xlim([0, None])
    ax.title.set_text("Ligand + Flexible Residues")
    ax.legend()

    plt.suptitle("Ligand RMSD Distribution")
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.pdf"))
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.png"))

    # rmsd = get_lig_rmsd_for_rank(1, df)
    # N = len(rmsd)
    # n = len(rmsd[rmsd < 1])
    # print(f"Number of sub-1Å top poses: {n} ({n / N * 100:.2f}%)")
