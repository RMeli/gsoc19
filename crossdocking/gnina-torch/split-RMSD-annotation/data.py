import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

df = pd.read_csv("../../cd-downsampled/analysis/SPLITflex1nc_annotated.csv")

sns.jointplot(x=df.rmsd, y=df.flexrmsd, kind="hex")
plt.savefig("plots/rmsd_vd_flexrmsd.png")
plt.savefig("plots/rmsd_vd_flexrmsd.pdf")

ligmax = 15
recmax = 4

print(
    df[(df.rmsd > ligmax) | (df.flexrmsd > recmax)].drop(
        columns=["score", "group", "lig_annotation", "rec_annotation", "obrmsd"]
    )
)

sns.jointplot(x=df[df.rmsd <= ligmax].rmsd, y=df[df.flexrmsd <= 5].flexrmsd, kind="hex")
plt.vlines(x=2, ymin=0, ymax=recmax, colors="grey", linestyles="dashed")
plt.hlines(y=1, xmin=0, xmax=ligmax, colors="grey", linestyles="dashed")
plt.savefig("plots/rmsd_vd_flexrmsd-zoom.png")
plt.savefig("plots/rmsd_vd_flexrmsd-zoom.pdf")

sns.jointplot(
    x=df[(df.rmsd <= ligmax) & (df["rank"] != 0)].rmsd,
    y=df[(df.flexrmsd <= recmax) & (df["rank"] != 0)].flexrmsd,
    kind="hex",
)
plt.vlines(x=2, ymin=0, ymax=recmax, colors="grey", linestyles="dashed")
plt.hlines(y=1, xmin=0, xmax=ligmax, colors="grey", linestyles="dashed")
plt.savefig("plots/rmsd_vd_flexrmsd-zoom-nocrystal.png")
plt.savefig("plots/rmsd_vd_flexrmsd-zoom-nocrystal.pdf")

# sns.heatmap(
#     data=df[["lig_annotation", "rec_annotation"]]
# )
# plt.savefig("plots/lig_rec_annotation.png")
# plt.savefig("plots/lig_rec_annotation.pdf")

g = sns.catplot(
    data=df[["lig_annotation", "rec_annotation"]],
    kind="count",
    x="rec_annotation",
    col="lig_annotation",
)
plt.tight_layout()
plt.savefig("plots/balance.png")
plt.savefig("plots/balance.pdf")

# rec_annotation for good ligand pose look very similar
# get actual counts
print(df[["rec_annotation"]][df["lig_annotation"] == 1].value_counts())
