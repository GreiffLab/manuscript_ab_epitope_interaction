import pandas as pd
import numpy as np

from util import compute_error, split_to_train_test, store_predictions, make_predictions_path


def compute_conditional_proba(train: pd.DataFrame, target: str, input: str) -> dict:
    cond_probs = train.groupby(target).agg(lambda x: pd.Series.mode(x)[0]).to_dict()[input]
    return cond_probs


def run(dataset: pd.DataFrame, split_count, training_percentage: float, target: str, input: str, result_path: str, comparison_fn) -> pd.DataFrame:
    performances = []
    for i in range(1, split_count + 1):
        train, test = split_to_train_test(dataset, training_percentage)
        mapping = compute_conditional_proba(train, target, input)
        performance = compute_error(test, mapping, target, input, comparison_fn)
        performances.append(performance)
        store_predictions(test, mapping, target, input, make_predictions_path(result_path, i))
    result_df = pd.DataFrame({"split_error": performances})
    result_df.to_csv(result_path)
    print("Conditional probabilities: mean error: {}, error std: {}".format(np.mean(performances), np.std(performances, dtype=np.float64)))
    return result_df
