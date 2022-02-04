import pandas as pd

rmsds = pd.read_csv("carlos_cd/rmsds2.csv")

missing = rmsds[rmsds["flexrmsd"].isna() | rmsds["fmaxrmsd"].isna()]
missing = missing[["pocket","protein","ligand"]]
missing = missing.drop_duplicates(ignore_index=True)

missing.to_csv("missing.csv", header=None, index=None)