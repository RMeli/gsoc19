import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

import argparse as ap

parser = ap.ArgumentParser(description="Compute center of geometry of a molecule.")
parser.add_argument("fnames", type=str, nargs="+", help="Full test results")
parser.add_argument("-o", "--output", type=str, default=None)
parser.add_argument("-g", "--graph", type=str, default=None)
args = parser.parse_args()

df_RMSD = pd.DataFrame(columns=["name", "RMSD", "Fold", "type", "method"])

for fold, fname in enumerate(args.fnames):

    df = pd.read_csv(fname)

    cnn_score = df.loc[df.groupby(["system"])["cnn_score"].idxmax()].drop(columns=["smina_score", "smina_rank"])
    smina_score = df.loc[df.groupby(["system"])["smina_score"].idxmin()].drop(columns=["cnn_score", "smina_rank"])
    lig_rmsd = df.loc[df.groupby(["system"])["lig_rmsd"].idxmin()].drop(columns=["cnn_score", "smina_score", "smina_rank", "fmax_rmsd"])
    fmax_rmsd = df.loc[df.groupby(["system"])["fmax_rmsd"].idxmin()].drop(columns=["cnn_score", "smina_score", "smina_rank", "lig_rmsd"])


    df_merged = cnn_score.merge(smina_score, on="system", suffixes=("_cnn","_smina"))
    df_merged = df_merged.merge(lig_rmsd, on="system")
    df_merged = df_merged.merge(fmax_rmsd, on="system")

    #print(df_merged.head())

    df_rmsd = df_merged.drop(columns=["cnn_score", "smina_score","system"]).melt(var_name="name", value_name="RMSD")
    df_rmsd["Fold"] = fold

    def t(row):
        if row["name"].split("_")[0] == "lig":
            return "Ligand"
        elif row["name"].split("_")[0] == "fmax":
            return "Flex"
        else:
            raise ValueError

    def m(row):
        if row["name"].split("_")[-1] == "smina":
            return "smina"
        elif row["name"].split("_")[-1] == "cnn":
            return "cnn"
        else:
            return "best"

    df_rmsd["type"] = df_rmsd.apply(t, axis=1)
    df_rmsd["method"] = df_rmsd.apply(m, axis=1)

    #print(df_rmsd.head())

    df_RMSD = pd.concat([df_RMSD, df_rmsd])

print(df_RMSD)

if args.graph == "boxen":
    plt.figure()
    g = sns.catplot(x="Fold", y="RMSD", hue="method", col="type", data=df_RMSD, hue_order=["smina", "cnn", "best"], kind="boxen", legend=False)
    g.set_titles("{col_name}")
    plt.legend()
    plt.suptitle("Top Pose RMSD Distributions")
    plt.tight_layout()
    plt.savefig(args.output)
elif args.graph == "distplot":
    plt.figure()
    bins = np.linspace(0, max(df_RMSD["RMSD"]), 150)
    g = sns.FacetGrid(df_RMSD.drop("name", axis=1), col="Fold", row="type", hue="method")
    g.map(sns.distplot, "RMSD", bins=bins)
    plt.legend()
    plt.xlim([0,6])
    plt.savefig(args.output)
else:
    raise Exception(f"{args.graph} is not a valid option.")