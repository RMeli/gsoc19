import pandas as pd

from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import OrdinalEncoder

def annotation(row):
    return 1 if row.rmsd <= 2.0 else 0

def ligname(row):
    pocket=row["pocket"]
    lig=row["lig"]
    rannk=row["rank"]
    
    return

def recname(row):
    return


df = pd.read_csv("analysis/rmsd_clean.csv", index_col=False)

# Assign different group to different pockets
# TODO: Assesss pocket similarity?
encoder = OrdinalEncoder()
df["group"] = encoder.fit_transform(df["pocket"].to_numpy().reshape(-1, 1))
df["group"] = df["group"].astype(int)

# Check that the number of unique pockets matches the number of groups
assert len(df["pocket"].unique()) == df["group"].max() + 1

# Split dataset in three folds (cross-validation)
# Randomly shuffle the data and keep the same pocket in the same fold
cv = GroupShuffleSplit(n_splits=3, random_state=42)
for fold, (train_idx, test_idx) in enumerate(
    cv.split(df.drop(columns=["group"]).to_numpy(), groups=df["group"].to_numpy())
):
    with open(f"files/train_{fold}.types", "w") as trout:
    for tr in train_idx:
        row =  df.iloc[tr]

        a = annotation(row)

        

