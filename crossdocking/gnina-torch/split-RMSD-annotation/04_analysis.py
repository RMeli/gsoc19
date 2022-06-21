# Plotly does not work on headless nodes!

import pandas as pd
import numpy as np
import plotly.express as px

df = pd.read_csv("all.csv")

df = df[(df.tt == "test") & (df.Epoch == 500)]
df = df.drop(columns=["Epoch", "tt"])
df = df.drop(columns=["Accuracy", "Loss (pose)", "Accuracy (flex)", "Loss (flex pose)"])
print(df.head())

hyperparameters = ["annotation", "model", "clustering", "stratification", "scale_flexpose_loss"]

mean = df.groupby(
    hyperparameters,
    as_index=False,
).mean()
assert np.allclose(mean.fold, np.ones_like(mean.fold))
mean = mean.drop(columns="fold")

print(mean.head())

# pd.plotting.parallel_coordinates(mean.drop(columns=["Balanced accuracy", "Balanced accuracy (flex)", "ROC AUC"]), "ROC AUC (flex)")
fig = px.parallel_categories(
    mean.drop(columns=["Balanced accuracy", "Balanced accuracy (flex)", "ROC AUC"]),
    color="ROC AUC (flex)",
    dimensions=hyperparameters,
)
# Scale parameter increases resolution
fig.write_image("pcord-flex.png", scale=5)

fig = px.parallel_categories(
    mean.drop(columns=["Balanced accuracy", "Balanced accuracy (flex)", "ROC AUC (flex)"]),
    color="ROC AUC",
    dimensions=hyperparameters,
)
# Scale parameter increases resolution
fig.write_image("pcord-lig.png", scale=5)