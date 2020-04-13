import os
import sys
from glob import glob

import pandas as pd
import numpy as np
from editdistance import eval as edit_distance

from marginal_proba import run as run_marginal_proba_prediction
from cond_proba import run as run_cond_proba_prediction
from cond_proba_with_prior import run as run_cond_proba_prior_prediction


def run_for_setting(dataset, target, input, comparison_fn, result_path, training_percentage, split_count, prior_strength):

    filename = "split_{}_training_{}_target_{}_comp_{}.csv".format(split_count, training_percentage, target, comparison_fn["name"])
    marginal_result = run_marginal_proba_prediction(dataset, split_count, training_percentage, target, input,
                                                    result_path + "marginal_proba_{}".format(filename),
                                                    comparison_fn["fn"])
    cond_proba_result = run_cond_proba_prediction(dataset, split_count, training_percentage, target, input,
                                                  result_path + "cond_proba_{}".format(filename),
                                                  comparison_fn["fn"])
    cond_proba_prior_result = run_cond_proba_prior_prediction(dataset, split_count, training_percentage, target, input,
                                                              result_path + "cond_proba_prior_{}".format(filename),
                                                              prior_strength=prior_strength, comparison_fn=comparison_fn["fn"])

    report = {"approach": ["marginal_proba", "cond_proba", "cond_proba_with_prior"],
              "error_mean": [marginal_result.mean().values[0], cond_proba_result.mean().values[0], cond_proba_prior_result.mean().values[0]],
              "error_standard_deviation": [marginal_result.std().values[0], cond_proba_result.std().values[0], cond_proba_prior_result.std().values[0]]}

    pd.DataFrame(report).to_csv("{}summary_{}".format(result_path, filename), index=False)


def prepare_dataset_paths(dataset_root_path: str = None):
    if dataset_root_path:
        dataset_paths = list(glob(dataset_root_path + "/*.tsv"))
    else:
        dataset_paths = list(glob("./dataset*/*.tsv"))
    return dataset_paths


def prepare_result_path(result_path, dataset_path, randomized_pairs: bool = False):
    path = result_path + str(os.path.basename(os.path.dirname(dataset_path))) + "/" \
           + str(os.path.basename(dataset_path).split('.')[0]) + "/"

    if randomized_pairs:
        path += "randomized_pairs/"

    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

    return path


def run(dataset_root_path: str = None):
    split_count = 10
    prediction_targets = ["epitope", "paratope"]
    training_percentage = 0.8
    result_path = "./results/"
    dataset_paths = prepare_dataset_paths(dataset_root_path)
    fn1 = {"name": "exact_match", "fn": lambda x, y: x != y}
    fn2 = {"name": "LD", "fn": lambda x, y: edit_distance(x, y) / max(len(x), len(y))}
    comparison_fns = [fn1, fn2]
    prior_strength = 0.001

    for dataset_path in dataset_paths:
        for target in prediction_targets:
            for comparison_fn in comparison_fns:

                # methods:
                dataset = pd.read_csv(dataset_path, sep="\t").dropna(axis=0)
                run_for_setting(dataset, target, prediction_targets[not prediction_targets.index(target)], comparison_fn,
                                prepare_result_path(result_path, dataset_path), training_percentage,
                                split_count, prior_strength)

                # randomized:
                dataset[target] = np.random.permutation(dataset[target].values)
                run_for_setting(dataset, target, prediction_targets[not prediction_targets.index(target)], comparison_fn,
                                prepare_result_path(result_path, dataset_path, True), training_percentage,
                                split_count, prior_strength)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run(sys.argv[1])
    else:
        run()
