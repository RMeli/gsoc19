import pandas as pd
import os
from collections import defaultdict

clusters = defaultdict(list)

cluster = 0
with open("new_pbis_3.5.sim", "r") as fin:
    for line in fin:
        # print(line)
        for pocket in line.strip().split():
            pocket = pocket.split("_")[0]
            # print(pocket, cluster)

            clusters["pocket"].append(pocket)
            clusters["cluster"].append(cluster)

        cluster += 1

# Check that all pockets have been assigned to a cluster
pockets = [
    d
    for d in os.listdir("../../carlos_cd")
    if os.path.isdir(os.path.join("../../carlos_cd", d))
]
# print(pockets)
assert len(pockets) == 92  # Check there are exactly 92 pockets
for p in pockets:
    if not p in clusters["pocket"]:
        print(f"Pocket {p} missing!")

df = pd.DataFrame(clusters)

carlosclusters = {}
for p in pockets:
    print(df[df.pocket == p].cluster)
