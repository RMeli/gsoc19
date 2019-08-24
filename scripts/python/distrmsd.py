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


def get_rmsd_for_rank(rank: int, rmsd_name: str, df: pd.DataFrame) -> np.ndarray:

    return df[rmsd_name].loc[df["rank"] == rank].values


def plot(df: pd.DataFrame, rmsd_name: str, maxrank: int, bins: int, ax):

    for rank in range(1, maxrank + 1):

        rmsd = get_rmsd_for_rank(rank, rmsd_name, df)

        sns.distplot(rmsd, bins=bins, hist=True, kde=True, label=f"rank {rank}", ax=ax)


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

    f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))

    plot(df, "rmsd_lig", args.maxrank, args.bins, ax1)
    ax1.set_xlabel("RMSD (Å)")
    ax1.set_xlim([0, 20])
    ax1.title.set_text("Ligand")
    ax1.axvline(x=2, color="k", linestyle="--", label="2 Å")
    ax1.legend()

    plot(df, "rmsd_flex", args.maxrank, args.bins, ax2)
    ax2.set_xlabel("RMSD (Å)")
    ax2.set_xlim([0, 3])
    ax2.title.set_text("Flexible Residues")
    ax2.axvline(x=1, color="k", linestyle="--", label="1 Å")
    ax2.legend()

    plot(df, "rmsd_fmax", args.maxrank, args.bins, ax3)
    ax3.set_xlabel("RMSD (Å)")
    ax3.set_xlim([0, 2.5])
    ax3.title.set_text("MAX Flexible Residue")
    ax3.axvline(x=1, color="k", linestyle="--", label="1 Å")
    ax3.legend()

    plt.suptitle("RMSD Distributions")
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.pdf"))
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.png"))

    ###

    rmsd_lig = get_rmsd_for_rank(1, "rmsd_lig", df)
    N = len(rmsd_lig)
    n = len(rmsd_lig[rmsd_lig < 2])
    print(f"Number of sub-2Å RMSD ligands (top pose): {n} ({n / N * 100:.2f}%)")

    rmsd_flex = get_rmsd_for_rank(1, "rmsd_flex", df)
    N = len(rmsd_flex)
    n = len(rmsd_flex[rmsd_flex < 1])
    print(
        f"Number of sub-1Å RMSD flexible residues (top pose): {n} ({n / N * 100:.2f}%)"
    )

    rmsd_fmax = get_rmsd_for_rank(1, "rmsd_fmax", df)
    N = len(rmsd_fmax)
    n = len(rmsd_fmax[rmsd_fmax < 1])
    print(
        f"Number of sub-1Å MAX RMSD of flexible residues (top pose): {n} ({n / N * 100:.2f}%)"
    )