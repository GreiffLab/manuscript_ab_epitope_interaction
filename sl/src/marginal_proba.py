import pandas as pd
import numpy as np

from util import split_to_train_test, store_predictions, make_predictions_path


def compute_max_marginal_proba(train: pd.DataFrame, predicting: str) -> str:
    max_marginal_proba = train[predicting].mode()[0]
    return max_marginal_proba


def assess_performance(df: pd.DataFrame, max_marginal_proba: str, target: str, fn):
    result = df.apply(lambda row: fn(max_marginal_proba, row[target]), axis=1, result_type="reduce")
    return sum(result)/len(result)


def run(dataset: pd.DataFrame, split_count, training_percentage: float, target: str, input: str, result_path: str, comparison_fn) -> pd.DataFrame:
    performances = []
    for i in range(1, split_count + 1):
        train, test = split_to_train_test(dataset, training_percentage)
        max_marginal_proba = compute_max_marginal_proba(train, target)
        performance = assess_performance(test, max_marginal_proba, target, comparison_fn)
        performances.append(performance)
        store_predictions(test, {"*": max_marginal_proba}, target, input, make_predictions_path(result_path, i))
    result_df = pd.DataFrame({"split_error": performances})
    result_df.to_csv(result_path)
    print("Marginal probabilities: mean error: {}, error std: {}".format(np.mean(performances), np.std(performances, dtype=np.float64)))
    return result_df

