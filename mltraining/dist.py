import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

import argparse as ap

parser = ap.ArgumentParser(description="Compute center of geometry of a molecule.")
parser.add_argument("fname", type=str, help="Full test results")
parser.add_argument("-o", "--output", type=str, default=None)
args = parser.parse_args()

df = pd.read_csv(args.fname)

cnn_score = df.loc[df.groupby(["system"])["cnn_score"].idxmax()].drop(columns=["smina_score", "smina_rank"])
smina_score = df.loc[df.groupby(["system"])["smina_score"].idxmin()].drop(columns=["cnn_score", "smina_rank"])
lig_rmsd = df.loc[df.groupby(["system"])["lig_rmsd"].idxmin()].drop(columns=["cnn_score", "smina_score", "smina_rank", "fmax_rmsd"])
fmax_rmsd = df.loc[df.groupby(["system"])["fmax_rmsd"].idxmin()].drop(columns=["cnn_score", "smina_score", "smina_rank", "lig_rmsd"])


df_merged = cnn_score.merge(smina_score, on="system", suffixes=("_cnn","_smina"))
df_merged = df_merged.merge(lig_rmsd, on="system")
df_merged = df_merged.merge(fmax_rmsd, on="system")

print(df_merged.head())

df_lig = df_merged[["lig_rmsd_cnn", "lig_rmsd_smina", "lig_rmsd"]].rename(columns={"lig_rmsd_cnn":"gnina", "lig_rmsd_smina" : "smina", "lig_rmsd" : "best"})
df_flex = df_merged[["fmax_rmsd_cnn", "fmax_rmsd_smina", "fmax_rmsd"]].rename(columns={"fmax_rmsd_cnn":"gnina", "fmax_rmsd_smina" : "smina", "fmax_rmsd" : "best"})

print(df_lig.head())
print(df_flex.head())

bin_size=0.1

plt.figure()
sns.distplot(df_lig["gnina"].to_numpy(), bins=np.arange(0, df_lig["gnina"].max() + bin_size, bin_size), label="CNN")
sns.distplot(df_lig["smina"].to_numpy(), bins=np.arange(0, df_lig["smina"].max() + bin_size, bin_size), label="smina")
sns.distplot(df_lig["best"].to_numpy(), bins=np.arange(0, df_lig["best"].max() + bin_size, bin_size), label="best")
plt.title("Top Pose RMSD Distributions for Ligand")
plt.xlabel("RMSD (A)")
plt.xlim([0,10])
plt.legend()
if args.output is None:
    plt.plot()
else:
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
    plt.savefig(args.output + "_flex.pdf")