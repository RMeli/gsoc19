import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt

prefix = "flex1"

dfs = []
for i in range(3):
    for stage in ["test", "train"]:
        df = pd.read_csv(f"training/{prefix}/training{i}_metrics_{stage}.csv")
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
    "Loss",
]

for tp in to_plot:
    df_to_plot = df[["Epoch", "Fold", "Phase", tp]]

    print(df_to_plot)

    plt.figure()
    sns.lineplot(
        x="Epoch", y=tp, hue="Phase", style="Phase", data=df_to_plot, markers=True
    )
    plt.legend(loc="lower center")
    plt.savefig(
        f"plots/{prefix}_{tp.replace('(','').replace(')','').replace(' ','_')}.png"
    )
