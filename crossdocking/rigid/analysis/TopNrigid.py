import pandas as pd

import sys
sys.path.insert(1, "../../cd-downsampled/analysis")

from TopN import process

if __name__ == "__main__":
    models = ["default2017", "default2018", "dense"]
    prefixes = ["cluster", "nc"]
    annotations = [""]

    df_score_list = []
    df_top_list = []
    for model in models:
        for prefix in prefixes:
            for annotation in annotations:
                # Not using the FLEX prefix would test the rigid docking CNN on the 
                # rifgid docking dataset. However the TopN script has an assertion on
                # the number of pockets that would fail because one pocket is discarded
                # from flexible docking for lack of actives
                df_score, df_top = process(model, annotation, prefix, prediction_prefix="FLEX")

                df_score_list.append(df_score)
                df_top_list.append(df_top)

    pd.concat(df_score_list).to_csv("allscores.csv", float_format="%.5f", index=False)
    df_top  = pd.concat(df_top_list)
    df_top.to_csv("allTopN.csv", float_format="%.5f", index=False)