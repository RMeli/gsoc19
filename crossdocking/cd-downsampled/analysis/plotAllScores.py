import seaborn as sns
import pandas as pd

from matplotlib import pyplot as plt

df = pd.read_csv("allscores.csv")

sns.catplot(data=df, x="label", col="annotation", kind="count", col_wrap=2)
plt.tight_layout()
plt.savefig("count.png")