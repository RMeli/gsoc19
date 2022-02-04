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
        "-b", "--bin_size", type=float, default=None, help="Histogram bins size"
    )
    parser.add_argument(
        "-opath", "--outputpath", type=str, default="", help="Output files path"
    )

    parser.add_argument(
        "--ligrmsd", type=str, default="rmsd_lig", help="Column name with ligand RMSD"
    )
    parser.add_argument(
        "--flexrmsd",
        type=str,
        default="rmsd_flex",
        help="Column name with side-chains RMSD",
    )
    parser.add_argument(
        "--flexmaxrmsd",
        type=str,
        default="rmsd_fmax",
        help="Column name with maximum side-chain RMSD",
    )

    return parser.parse_args(args)


def get_rmsd_for_rank(rank: int, rmsd_name: str, df: pd.DataFrame) -> np.ndarray:

    return df[rmsd_name].loc[df["rank"] == rank].values


def plot(df: pd.DataFrame, rmsd_name: str, maxrank: int, bin_size: int, ax):
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

    for rank in range(0, maxrank + 1):

        rmsd = get_rmsd_for_rank(rank, rmsd_name, df)

        l = f"smina rank {rank}" if rank > 0 else "smina crystal"

        bins = np.arange(rmsd.min(), rmsd.max(), bin_size)

        sns.histplot(rmsd, bins=bins, kde=True, label=l, color=colors[rank], ax=ax)


if __name__ == "__main__":

    import os

    args = parse()

    df = pd.read_csv(args.input)

    # Remove swap flexrmsd with flexobrmsd column, if present
    # For some systems, spyrmsd is too slow for flexible residues
    try:
        df = df.drop(columns=["flexrmsd"])
        df = df.rename(columns={"flexobrmsd": "flexrmsd"})
    except KeyError:
        # No flexobrmsd value
        # In the original GSoC project there was no spyrmsd value
        pass

    bad_values = [np.nan, np.inf, -np.inf]
    if df.isin(bad_values).sum(axis=0).max() > 0:
        print(f"WARNING: {args.input} contains NaN or Inf:")
        print(df.isin(bad_values).sum(axis=0).max())
        bidx_lig = df.index[
            df[args.ligrmsd].isin(bad_values)
        ]  # Get row indices of bad values

        # Check there are no bad values for the ligand
        assert len(df.iloc[bidx_lig]) == 0

        bidx_flex = df.index[
            df[args.flexrmsd].isin(bad_values)
        ]  # Get row indices of bad values
        print(df.iloc[bidx_flex])

        # Replace infinite values with NaNs and drop them
        df = df.replace([np.inf, -np.inf], np.nan)  # Change Inf with NaN
        df = df.dropna()  # Drop rows with NaN

    f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))

    plot(df, args.ligrmsd, args.maxrank, args.bin_size, ax1)
    ax1.set_xlabel("RMSD (Å)")
    ax1.set_ylabel("Normalised density")
    ax1.set_xlim([0, 15])
    ax1.title.set_text("Ligand")
    ax1.axvline(x=2, color="k", linestyle="--", label="2 Å")
    ax1.legend()

    plot(df, args.flexrmsd, args.maxrank, args.bin_size, ax2)
    ax2.set_xlabel("RMSD (Å)")
    ax2.set_xlim([0, 3])
    ax2.title.set_text("Flexible Residues")
    handle = ax2.axvline(x=1, color="k", linestyle="--", label="1 Å")
    ax2.legend()
    # ax2.legend(handles=[handle])

    plot(df, args.flexmaxrmsd, args.maxrank, args.bin_size, ax3)
    ax3.set_xlabel("RMSD (Å)")
    ax3.set_xlim([0, 4])
    ax3.title.set_text("MAX Flexible Residue")
    handle = ax3.axvline(x=1, color="k", linestyle="--", label="1 Å")
    ax3.legend()
    # ax3.legend(handles=[handle])

    # plt.suptitle("RMSD Distributions")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.pdf"))
    plt.savefig(os.path.join(args.outputpath, f"distrmsd.png"))

    rmsd_lig = get_rmsd_for_rank(1, args.ligrmsd, df)
    N = len(rmsd_lig)
    n = len(rmsd_lig[rmsd_lig < 2])
    print(f"Number of sub-2Å RMSD ligands (top pose): {n} ({n / N * 100:.2f}%)")

    rmsd_flex = get_rmsd_for_rank(1, args.flexrmsd, df)
    N = len(rmsd_flex)
    n = len(rmsd_flex[rmsd_flex < 1])
    print(
        f"Number of sub-1Å RMSD flexible residues (top pose): {n} ({n / N * 100:.2f}%)"
    )

    rmsd_fmax = get_rmsd_for_rank(1, args.flexmaxrmsd, df)
    N = len(rmsd_fmax)
    n = len(rmsd_fmax[rmsd_fmax < 1])
    print(
        f"Number of sub-1Å MAX RMSD of flexible residues (top pose): {n} ({n / N * 100:.2f}%)"
    )
