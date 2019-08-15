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
        "-opath", "--outputpath", type=str, default="", help="Output files path"
    )

    return parser.parse_args(args)


def get_num_for_threshold(df: pd.DataFrame, threshold: float, rmsd_name: str) -> float:

    # Number of rows
    n = len(df.index)

    n_rmsd = df[rmsd_name].loc[df[rmsd_name] < threshold].count()

    return n_rmsd / n * 100.0


def plot(df: pd.DataFrame, rmsd_name: str, ax, tmin: float = 0, tmax: float = 4):

    threshold = np.linspace(tmin, tmax, 100)

    num = []
    for t in threshold:
        num.append(get_num_for_threshold(df, t, rmsd_name))

    num = np.array(num)

    sns.lineplot(threshold, num, ax=ax)


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

    plot(df, "rmsd_lig", ax1, 0, 10)
    ax1.set_xlabel("RMSD Threshold (Å)")
    ax1.set_ylabel("% of Good Poses")
    ax1.set_xlim([0, None])
    ax1.title.set_text("Ligand")

    plot(df, "rmsd_flex", ax2, 0, 3)
    ax2.set_xlabel("RMSD Threshold (Å)")
    ax2.set_xlim([0, None])
    ax2.title.set_text("Flexible Residues")

    plot(df, "rmsd_fmax", ax3, 0, 5)
    ax3.set_xlabel("RMSD Threshold (Å)")
    ax3.set_xlim([0, None])
    ax3.title.set_text("MAX Flexible Residue")

    plt.suptitle("Number of Good Poses vs RMSD Threshold")
    plt.savefig(os.path.join(args.outputpath, f"rmsdthreshold.pdf"))
    plt.savefig(os.path.join(args.outputpath, f"rmsdthreshold.png"))

    plt.show()
