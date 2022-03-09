import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Plot results of training.")
parser.add_argument(
    "model", type=str, choices=["default2017-nc", "default2018-nc", "dense-nc"]
)
parser.add_argument("prefix", type=str, choices=["flex1", "felex2", "max2"])

args = parser.parse_args()

dfs = []
for i in range(3):
    for stage in ["test", "train"]:
        df = pd.read_csv(
            f"training/{args.prefix}/{args.model}/training{i}_metrics_{stage}.csv"
        )
        df["Fold"] = i
        df["Phase"] = stage

        dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

to_plot = [
    "ROC AUC",
    "ROC AUC (flex)",
    "Balanced accuracy",
    "Balanced accuracy (flex)",
    "Loss (pose)",
    "Loss (flex pose)",
]

for tp in to_plot:
    df_to_plot = df[["Epoch", "Fold", "Phase", tp]]

    print(df_to_plot)

    plt.figure()
    sns.lineplot(
        x="Epoch", y=tp, hue="Phase", style="Phase", data=df_to_plot, markers=True
    )
    plt.legend(loc="lower center")

    figname = f"plots/{args.prefix}_{args.model}_{tp.replace('(','').replace(')','').replace(' ','_')}"
    plt.savefig(
        f"{figname}.png"
    )
    plt.savefig(
        f"{figname}.pdf"
    )
