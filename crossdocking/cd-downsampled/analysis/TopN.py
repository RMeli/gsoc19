import pandas as pd
import numpy as np


def load_scores(df, model: str, annotation: str, prefix: str, prediction_prefix=""):
    """
    Load CNN score from prediction files and create a clean dataframe together with
    RMSD data.

    Parameters
    ----------
    df: pd.DataFrame
        Data frame with RMSD values
    model: str
        Name of the model
    annotation: str
        Annotation used to train the CNN
    prefix: str
        Prefix to distinguish between clustered and not clustered
    pprediction_predix:
        Optional prediction prefix for rigid docking CNN applied to flexible dataset

    Returns
    -------
    pd.DataFrame
        Dataframe with full information, ready for analysis

    Notes
    -----
    The CNN scoring function trained with systems from rigid docking can be used to
    to score poses obtained with flexible docking. When this is done, the FLEX prefix
    is appended to the output prediction files.
    """

    tfname = annotation if annotation else "ligand"

    print(tfname, model, prefix)

    # Load predicted scores
    dfs = []
    for i in range(3):
        dfs.append(
            pd.read_csv(
                f"../training/{tfname}/{model}/{prediction_prefix}{annotation}{prefix}test{i}.out",
                sep=" ",
                header=None,
            )
        )

    score = pd.concat(dfs)

    if prefix == "nc":  # No affinity column
        to_drop = [1, 2, 3, 5, 6, 7]
    elif prefix == "cluster":
        to_drop = [1, 2, 3, 4, 6, 7, 8]

    # Drop useless columns
    score.drop(columns=to_drop, inplace=True)

    # Rename either column 4 (no clustering) or 5 (clustering) as ligand name
    # The shift between the two cases is given by the presence/absence of the affiniy
    # column, used here for annotating different clusters
    score.rename(columns={0: "CNNscore", 4: "ligname", 5: "ligname"}, inplace=True)

    score.dropna(inplace=True)  # Last row contains NaN (it is  actually a comment)

    # CNNscore is float since the last row contained # in the column (comment)
    score["CNNscore"] = score["CNNscore"].astype(float)

    # Exrtract information from ligand name (pocket, protein, ligand, rank)
    try:
        # Flexible docking is based on the carlos_cd dataset of the following form
        # DATA_FOLDER/POCKET/PDB_Structures/XXXX_PRO_YYYY_LIG_*_pRANK.gninatypes
        score[["pocket", "ligname"]] = score["ligname"].str.split("/", expand=True)[
            [1, 3]
        ]
    except KeyError:
        # Rigid docking is based on the wierbowski dataset
        # The dataset is the same as carlos_cd, but renamed and without the PDB_Structures
        # DATA_FOLDER/POCKET/XXXX_PRO_YYYY_LIG_*_pRANK.gninatypes
        # This means that splitting the file path results in one less element
        # which corresponded to the PDB_Structures folder, now gone
        score[["pocket", "ligname"]] = score["ligname"].str.split("/", expand=True)[
            [1, 2]
        ]

    score[["protein", "ligand", "rank"]] = score["ligname"].str.split("_", expand=True)[
        [0, 2, 10]
    ]
    score["rank"] = (
        score["rank"].str.replace("p", "").str.replace(".gninatyes", "").astype(int)
    )
    score.drop(columns="ligname", inplace=True)

    df_score = df.merge(score, on=["pocket", "protein", "ligand", "rank"])

    # Store additional information about the run
    # This will allow to combine all data in a single DataFrame for easier analysis
    df_score["model"] = model.split("-")[0]
    df_score["annotation"] = tfname
    df_score["prefix"] = prefix

    return df_score


def topN(df, nmax: int):
    """
    Compute TopN metric (for N in [1, NMAX])

    Parameters
    ----------
    df: pd.DataFrame
        Data frame with RMSD and pose scores
    nmax: int
        Maximum value of N

    Returns
    -------
    np.ndarray
        Array with results

    Notes
    -----
    The TopN metric is the percentage of targets with a good pose amongst the top N
    ranked poses. In the case of cross-docking this metric is averaged by pocket.
    """

    n_pockets = 0

    # Store top N metrics
    top_smina = [0] * nmax
    top_gnina = [0] * nmax
    top_best = [0] * nmax

    # Loop over pockets
    for _, group in df.groupby(by=["pocket"]):
        n_smina = [0] * nmax
        n_gnina = [0] * nmax
        n_best = [0] * nmax

        # Compute percentage of targets with good pose in top N
        n_targets = 0
        for _, tgroup in group.groupby("protein"):
            smina = tgroup.sort_values(by="score", ascending=True)
            gnina = tgroup.sort_values(by="CNNscore", ascending=False)
            best = tgroup.sort_values(by="rmsd", ascending=True)

            for n in range(1, nmax + 1):
                # At least one good pose amongst the top N
                if (smina["label"].iloc[:n] == 1).any():
                    n_smina[n - 1] += 1

                if (gnina["label"].iloc[:n] == 1).any():
                    n_gnina[n - 1] += 1

                if (best["label"].iloc[:n] == 1).any():
                    n_best[n - 1] += 1

            n_targets += 1

        # Accumulate results for all targets
        for n in range(1, nmax + 1):
            top_smina[n - 1] += n_smina[n - 1] / n_targets * 100
            top_gnina[n - 1] += n_gnina[n - 1] / n_targets * 100
            top_best[n - 1] += n_best[n - 1] / n_targets * 100

        n_pockets += 1

    # One pocket has been removed for lack of actives
    assert n_pockets == 91 or n_pockets == 92

    # Return TopN of targets, averaged per pocket
    top_smina_avg = np.array(top_smina) / n_pockets
    top_gnina_avg = np.array(top_gnina) / n_pockets
    top_best_avg = np.array(top_best) / n_pockets

    t = np.array(
        [list(range(1, nmax + 1)), top_smina_avg, top_gnina_avg, top_best_avg]
    ).T

    return pd.DataFrame(t, columns=["N", "smina", "gnina", "best"])


def process(model, annotation, prefix, nmax=10, prediction_prefix=""):
    if prefix == "nc":
        modelname = f"{model}-noaffinity-nostratified"
    elif prefix == "cluster":
        modelname = f"{model}-noaffinity"

    df = pd.read_csv(f"{annotation}nc_annotated.csv")

    # !!!!!!!!!!
    # Change the column name "annotation" to "label"
    # The name annotation is also used for the flexible docking annotation
    # Label is more appropriate here anyways
    df.rename(columns={"annotation": "label"}, inplace=True)
    # !!!!!!!!!!

    df_score = load_scores(df, modelname, annotation, prefix, prediction_prefix)

    df_top_list = []

    # Compute TopN metric both with and without the crystal structure
    for crystal in ["yes", "no"]:
        if crystal == "yes":
            df_top = topN(df_score, nmax)
        elif crystal == "no":
            # Remove crystal
            df_top = topN(df_score[df_score["rank"] != 0], nmax)
        else:
            raise Exception

        df_top["annotation"] = annotation if annotation else "ligand"
        df_top["prefix"] = prefix
        df_top["model"] = model
        df_top["crystal"] = crystal

        df_top["N"] = df_top["N"].astype(int)

        df_top_list.append(df_top)

        # df_top.to_csv(f"TopN/{crystal}{modelname}_{annotation}{prefix}.csv", index=None)

    return df_score, pd.concat(df_top_list)


if __name__ == "__main__":
    models = ["default2017", "default2018", "dense"]
    # models = ["default2017", "default2018"]
    prefixes = ["cluster", "nc"]
    # prefixes = ["nc"]
    annotations = ["", "flex", "max1", "max2"]
    # annotations = ["",]

    df_score_list = []
    df_top_list = []
    for model in models:
        for prefix in prefixes:
            for annotation in annotations:
                df_score, df_top = process(model, annotation, prefix)

                df_score_list.append(df_score)
                df_top_list.append(df_top)

    pd.concat(df_score_list).to_csv("allscores.csv", float_format="%.5f", index=False)
    df_top = pd.concat(df_top_list)
    df_top.to_csv("allTopN.csv", float_format="%.5f", index=False)
