import pandas as pd
import seaborn as sns
import os

from matplotlib import pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Plot results of training.")
parser.add_argument(
    "output",
    type=str,
)
parser.add_argument("plot", type=str)
args = parser.parse_args()

wdir = os.path.dirname(args.output)

dfs = []
for i in range(3):
    for stage in ["test", "train"]:
        df = pd.read_csv(os.path.join(wdir, f"training{i}_metrics_{stage}.csv"))
        df["Fold"] = i
        df["Phase"] = stage

        dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

to_plot = {
    "roc-auc": "ROC AUC",
    "roc-auc-flex": "ROC AUC (flex)",
    "balanced-accuracy": "Balanced accuracy",
    "balanced-accuracy-flex": "Balanced accuracy (flex)",
    "loss-pose": "Loss (pose)",
    "loss-flex-pose": "Loss (flex pose)",
}

df_to_plot = df[["Epoch", "Fold", "Phase", to_plot[args.plot]]]

plt.figure()
sns.lineplot(
    x="Epoch",
    y=to_plot[args.plot],
    hue="Phase",
    style="Phase",
    data=df_to_plot,
    markers=True,
)
plt.legend(loc="lower center")

plt.savefig(args.output)
plt.savefig(args.output.replace(".png", ".pdf"))
