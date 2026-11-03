from __future__ import annotations

import pandas as pd

from app.evidence_intelligence import disagreement_probability, expected_information_gain, select_robust_method


def test_robust_selection_minimizes_worst_scheme_rank() -> None:
    methods = ("median_potency", "student_t_regression", "xgboost", "phytotrust_m0", "m1_no_protocol", "phytotrust_m1")
    metrics = {}
    for scheme, order in {
        "study_grouped": ("student_t_regression", "phytotrust_m0", "xgboost", "median_potency", "m1_no_protocol", "phytotrust_m1"),
        "forward_temporal": ("xgboost", "student_t_regression", "phytotrust_m0", "median_potency", "m1_no_protocol", "phytotrust_m1"),
    }.items():
        metrics[scheme] = {"all": {
            method: {"nll": float(order.index(method)), "mean_absolute_calibration_error": 0.1}
            for method in methods
        }}
    assert select_robust_method(metrics)["selected_method"] == "student_t_regression"


def test_disagreement_probability_rises_with_interval_width() -> None:
    narrow = pd.Series({"prediction": 6.0, "prior_median_pactivity": 6.0, "lower_90": 5.8, "upper_90": 6.2})
    wide = pd.Series({"prediction": 6.0, "prior_median_pactivity": 6.0, "lower_90": 4.0, "upper_90": 8.0})
    assert disagreement_probability(wide) > disagreement_probability(narrow)


def test_information_gain_increases_with_uncertainty_and_decreases_with_support() -> None:
    assert expected_information_gain(2.0, 1) > expected_information_gain(1.0, 1)
    assert expected_information_gain(2.0, 1) > expected_information_gain(2.0, 4)
