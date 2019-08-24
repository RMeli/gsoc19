import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

import argparse as ap

parser = ap.ArgumentParser(description="Compute center of geometry of a molecule.")
parser.add_argument("fnames", type=str, nargs="+", help="Full test results")
parser.add_argument("-o", "--output", type=str, default=None)
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

    print(df_merged.head())

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

    print(df_rmsd.head())

    df_RMSD = pd.concat([df_RMSD, df_rmsd])

g = sns.catplot(x="Fold", y="RMSD", hue="method", col="type", data=df_RMSD, kind="boxen", legend=False, aspect=0.5)
g.set_titles("{col_name}")
plt.legend()
plt.suptitle("Top Pose RMSD Distributions")
plt.tight_layout()
plt.savefig("dist.pdf")
plt.show()

#bin_size=0.1

#plt.figure()
#sns.distplot(df_lig["gnina"].to_numpy(), bins=np.arange(0, df_lig["gnina"].max() + bin_size, bin_size), label="CNN")
#sns.distplot(df_lig["smina"].to_numpy(), bins=np.arange(0, df_lig["smina"].max() + bin_size, bin_size), label="smina")
#sns.distplot(df_lig["best"].to_numpy(), bins=np.arange(0, df_lig["best"].max() + bin_size, bin_size), label="best")
#plt.title("Top Pose RMSD Distributions for Ligand")
#plt.xlabel("RMSD (A)")
#plt.xlim([0,10])
#plt.legend()
#if args.output is None:
#    plt.plot()
""" else:
    plt.savefig(args.output + "_lig.png")
    plt.savefig(args.output + "_lig.pdf")

plt.figure()
sns.distplot(df_flex["gnina"].to_numpy(), bins=np.arange(0, df_flex["gnina"].max() + bin_size, bin_size), label="CNN")
sns.distplot(df_flex["smina"].to_numpy(), bins=np.arange(0, df_flex["smina"].max() + bin_size, bin_size), label="smina")
sns.distplot(df_flex["best"].to_numpy(), bins=np.arange(0, df_flex["best"].max() + bin_size, bin_size), label="best")
plt.title("Top Pose RMSD Distributions for Flex")
plt.xlabel("RMSD (A)")
plt.xlim([0,6])
plt.legend()
if args.output is None:
    plt.plot()
else:
    plt.savefig(args.output + "_flex.png")
    plt.savefig(args.output + "_flex.pdf") """