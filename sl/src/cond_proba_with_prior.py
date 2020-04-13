import pandas as pd
import numpy as np

from util import split_to_train_test, compute_error, store_predictions, make_predictions_path


def compute_non_normalized_proba(counts, target_counts, target_index, prior_strength):
    non_normalized_proba = {}
    for combination_name in counts:
        non_normalized_proba[combination_name] = counts[combination_name] + prior_strength * target_counts[combination_name[target_index]]
    return non_normalized_proba


def compute_normalized_proba(counts, non_normalized_proba, input_index):
    normalized_proba = {}
    for combination_name in counts:
        normalized_proba[combination_name] = non_normalized_proba[combination_name] / \
                                        sum([non_normalized_proba[key] for key in non_normalized_proba
                                             if key[input_index] == combination_name[input_index]])
    return normalized_proba


def compute_max_proba(input_values, target_index, input_index, normalized_proba):
    max_probabilities = {}
    for input_value in input_values:
        prob_classes = {key[target_index]: normalized_proba[key] for key in normalized_proba if key[input_index] == input_value}
        assert np.isclose(sum(prob_classes[key] for key in prob_classes), 1)
        max_probabilities[input_value] = max(prob_classes, key=lambda key: prob_classes[key])
    return max_probabilities


def compute_probabilities(train: pd.DataFrame, target: str, input: str, prior_strength: float) -> dict:

    input_values = train[input].unique()

    counts = train.groupby(train.columns.to_list(), as_index=False).size().to_dict()
    target_index = train.columns.to_list().index(target)
    input_index = train.columns.to_list().index(input)
    target_counts = train.groupby(target).size().to_dict()

    # calculate non-normalized probabilities: tmp = count(X,y) + k * count(y)
    non_normalized_proba = compute_non_normalized_proba(counts, target_counts, target_index, prior_strength)

    # normalize probabilities: P(Y|X) = tmp / sum(tmp | X)
    normalized_proba = compute_normalized_proba(counts, non_normalized_proba, input_index)

    # construct max probability mapping: X: Y | Y = argmax Y (P(Y|X)|X)
    max_probabilities = compute_max_proba(input_values, target_index, input_index, normalized_proba)

    return max_probabilities


def run(dataset: pd.DataFrame, split_count: int, training_percentage: float, target: str, input: str, result_path: str, comparison_fn,
        prior_strength: float) -> pd.DataFrame:
    performances = []
    for i in range(1, split_count + 1):
        train, test = split_to_train_test(dataset, training_percentage)
        probabilities = compute_probabilities(train, target, input, prior_strength)
        performance = compute_error(test, probabilities, target, input, comparison_fn)
        performances.append(performance)
        store_predictions(test, probabilities, target, input, make_predictions_path(result_path, i))
    result_df = pd.DataFrame({"split_error": performances})
    result_df.to_csv(result_path)
    print("Conditional probabilities with prior: mean error: {}, error std: {}".format(np.mean(performances), np.std(performances, dtype=np.float64)))
    return result_df
