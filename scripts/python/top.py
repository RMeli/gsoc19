import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

def smina_results(df_results: pd.DataFrame, lig_only: bool = False, csv_name: str = "fulltest", n_folds: int=3):

    for fold in range(n_folds):
        # Read and cleanup fold results
        df = pd.read_csv(f"{csv_name}{fold}.csv")
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Get top smina score (lowest value)
        df_score = df.loc[df.groupby(["system"])["smina_score"].idxmin()]
        N = len(df_score)

        # Good poses as top
        if lig_only:
            df_good = df_score.loc[(df_score["lig_rmsd"] < 2.0)]
        else:
            df_good = df_score.loc[(df_score["lig_rmsd"] < 2.0) & (df_score["fmax_rmsd"] < 1)]
        n = len(df_good)

        # Percentage of good poses scored as top
        r = n / N * 100

        #print(f"SMINA - Top 1: {r:.5f}")

        df_results = df_results.append({"sf" : "smina", "top": r}, ignore_index=True)

    return df_results

def cnn_results(df_results: pd.DataFrame, lig_only: bool = False, csv_name: str = "fulltest", n_folds: int=3):

    for fold in range(n_folds):
        # Read and cleanup fold results
        df = pd.read_csv(f"{csv_name}{fold}.csv")
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Get top cnn score (highest value)
        df_score = df.loc[df.groupby(["system"])["cnn_score"].idxmax()]
        N = len(df_score)

        # Good poses as top
        if lig_only:
            df_good = df_score.loc[(df_score["lig_rmsd"] < 2.0)]

        else:
            df_good = df_score.loc[(df_score["lig_rmsd"] < 2.0) & (df_score["fmax_rmsd"] < 1)]
        n = len(df_good)

        # Percentage of good poses as top
        r = n / N * 100

        #print(f"CNN - Top 1: {r:.5f}")

        df_results = df_results.append({"sf" : "cnn", "top": r}, ignore_index=True)

    return df_results

def  best_results(df_results: pd.DataFrame, lig_only: bool = False, csv_name: str = "fulltest", n_folds: int=3):

    for fold in range(n_folds):
        # Read and cleanup fold results
        df = pd.read_csv(f"{csv_name}{fold}.csv")
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Good poses annotation
        if lig_only:
            df["good"] = (df["lig_rmsd"] < 2.0)
        else:
            df["good"] = (df["lig_rmsd"] < 2.0) & (df["fmax_rmsd"] < 1)

        # Count the number of good poses
        df_good = df.groupby(["system"])["good"].sum()

        n = len(df_good[df_good > 0.5])
        N = len(df_good)

        # Percentage of good poses as top
        r = n / N * 100

        #print(f"CNN - Top 1: {r:.5f}")

        df_results = df_results.append({"sf" : "best", "top": r}, ignore_index=True)

    return df_results

if __name__ == "__main__":

    import argparse as ap
    from typing import Optional

    def parse(args: Optional[str] = None) -> ap.Namespace:

        parser = ap.ArgumentParser(
            description="Good top pose."
        )

        parser.add_argument("csv_name", type=str, help="CSV result file name (without extension)")
        parser.add_argument("-n", "--n_folds", type=int, default=3, help="Number of folds")
        parser.add_argument("-o", "--output", type=str, default="top.pdf", help="Output graph")
        parser.add_argument("-l", "--ligonly", default=False, action="store_true", help="Ligand-based annotation")

        return parser.parse_args(args)

    args = parse()
    
    # SF: Scoring Function
    df_results = pd.DataFrame(columns=["sf", "top"])

    df_results = smina_results(df_results, args.ligonly, args.csv_name, args.n_folds)
    df_results = cnn_results(df_results, args.ligonly, args.csv_name, args.n_folds)
    df_results = best_results(df_results, args.ligonly, args.csv_name, args.n_folds)

    sns.boxplot(data=df_results, x="sf", y="top")
    sns.swarmplot(data=df_results, x="sf", y="top")
    plt.xlabel("")
    plt.ylabel("Percentage of Targets with GOOD Top Pose")
    plt.savefig(args.output)

