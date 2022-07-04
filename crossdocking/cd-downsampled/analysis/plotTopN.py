import seaborn as sns
import pandas as pd

from matplotlib import pyplot as plt

df = pd.read_csv("allTopN.csv")

df_melted = df.melt(
    value_vars=["best", "smina", "gnina"],
    id_vars=["N", "annotation", "prefix", "model", "crystal"],
)

print(df_melted)

# Rename GNINA column with the actual model name and rename columns for plotting
idxs = df_melted["variable"] == "gnina"
df_melted.loc[idxs, "variable"] = df_melted.loc[idxs, "model"]
df_melted.drop(columns="model", inplace=True)
df_melted.rename(columns={"variable": "method", "value": "TopN (%)"}, inplace=True)

print(df_melted)

# Rename prefix
df_melted["prefix"] = df_melted["prefix"].str.replace("cluster", "clustering")
df_melted["prefix"] = df_melted["prefix"].str.replace("nc", "no clustering")


def plot(df, name):
    g = sns.relplot(
        data=df,
        x="N",
        y="TopN (%)",
        hue="method",
        style="crystal",
        col="prefix",
        row="annotation",
        kind="line",
        markers=True,
    )
    g.tight_layout()  # Excludes legend; plt.tight_layout() includes legend
    plt.ylim([0, 100])
    plt.savefig(f"plots/TopN-{name}.png")
    plt.savefig(f"plots/TopN-{name}.pdf")


plot(df_melted[df_melted.annotation != "ligand"], "flexannotation")
plot(df_melted[df_melted.annotation == "ligand"], "ligannotation")
