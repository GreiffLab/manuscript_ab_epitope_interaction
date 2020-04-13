import os

import numpy as np
import pandas as pd


def split_to_train_test(df: pd.DataFrame, training_percentage: float):
    mask = np.random.randn(len(df)) < training_percentage
    train = df[mask]
    test = df[~mask]
    return train, test


def compute_error(df: pd.DataFrame, mapping: dict, target: str, input: str, fn):
    result = df.apply(lambda row: fn(mapping[row[input]], row[target])
                                  if row[input] in mapping else True, axis=1, result_type="reduce")
    return sum(result)/len(result)


def get_prediction(mapping: dict, row, input):
    if row[input] in mapping:
        return mapping[row[input]]
    elif "*" in mapping:
        return mapping["*"]
    else:
        return None


def store_predictions(df: pd.DataFrame, mapping: dict, target: str, input: str, filepath: str):

    predictions = df.apply(lambda row: pd.Series([row[input],
                                                  get_prediction(mapping, row, input),
                                                  row[target]],
                                                 index=["input_{}".format(input),
                                                        "prediction_{}".format(target),
                                                        "true_output_{}".format(target)]),
                           axis=1)

    predictions.to_csv(filepath, index=False)


def make_predictions_path(split_result_path: str, split: int):
    directory = os.path.dirname(split_result_path) + "/predictions/"

    if not os.path.isdir(directory):
        os.mkdir(directory)

    path = directory + str(os.path.basename(split_result_path)).split("_comp_")[0] + "_split_{}.csv".format(split)
    return path
