"""
Collect all CSV files from the different running folders in order to perform analysis
and visualization of all different runs, in order to explore the hyperparameter space.

This script produces a CSV file "all.csv" containing all the data.
"""

import pandas as pd
import os

folds = list(range(3))
models = ["default2017", "default2018", "dense"]
clusterings = ["nc", "cluster"]
annotations = ["flex05", "flex1", "flex2", "max2"]
stratify = ["0,1,0.5", "0,0,0"]
scale_flexpose_loss = [1.0, 5.0, 10.0]

traindir = "train"

# Collects all datastes together
dfs = []

for a in annotations:
    for m in models:
        for c in clusterings:
            for s in stratify:
                for l in scale_flexpose_loss:
                    for f in folds:
                        for t in ["train", "test"]:
                            fname = os.path.join(traindir, a, f"{m}-{c}", f"s{s}-floss{l}", f"training{f}_metrics_{t}.csv")
                            
                            df = pd.read_csv(fname)

                            df["annotation"] = a
                            df["model"] = m
                            df["clustering"] = c
                            df["stratification"] = s
                            df["scale_flexpose_loss"] = l
                            df["fold"] = f
                            df["tt"] = t

                            dfs.append(df)

df_all = pd.concat(dfs)
df_all.to_csv("all.csv", index=False)