from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from app.phytotrust_m0 import (
    HierarchicalStudentTM0,
    M0Config,
    _assert_fold_safe,
    cluster_bootstrap_differences,
    evaluate_predictions,
)


def synthetic_frame(seed: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for target_index, target in enumerate(("T1", "T2")):
        for study_index in range(12):
            shift = rng.normal(0, 0.45)
            for compound_index in range(3):
                prior = 5.5 + rng.normal(0, 0.5)
                outcome = prior + (0.2 if target_index else -0.1) + shift + rng.standard_t(5) * 0.18
                rows.append({
                    "target_chembl_id": target,
                    "heldout_study_key": f"{target}-S{study_index}",
                    "molecule_id": f"{target}-S{study_index}-C{compound_index}",
                    "heldout_year": 2000 + study_index,
                    "prior_median_pactivity": prior,
                    "heldout_median_pactivity": outcome,
                    "is_verified_natural_product": compound_index == 0,
                })
    return pd.DataFrame(rows)


def test_hierarchical_estimation_recovers_nonzero_study_variation() -> None:
    model = HierarchicalStudentTM0(M0Config(quadrature_nodes=11, optimizer_maxiter=400)).fit(synthetic_frame())
    parameters = model.parameters()
    assert parameters["residual_scale"] > 0
    assert all(value > 0.10 for value in parameters["target_between_study_sd"].values())
    assert all(value < 1.50 for value in parameters["target_between_study_sd"].values())


def test_unseen_study_uncertainty_is_integrated_not_conditioned() -> None:
    frame = synthetic_frame()
    model = HierarchicalStudentTM0(M0Config(quadrature_nodes=11, optimizer_maxiter=400)).fit(frame)
    first = frame.iloc[[0]].copy()
    altered = first.copy()
    altered["heldout_median_pactivity"] = 99.0
    prediction_a = model.predict_distribution(first)
    prediction_b = model.predict_distribution(altered)
    assert prediction_a["prediction"][0] == pytest.approx(prediction_b["prediction"][0])
    assert prediction_a["lower_80"][0] == pytest.approx(prediction_b["lower_80"][0])
    assert prediction_a["upper_80"][0] == pytest.approx(prediction_b["upper_80"][0])
    assert prediction_a["between_study_sd"][0] > 0


def test_hierarchical_intervals_include_between_study_uncertainty() -> None:
    frame = synthetic_frame()
    config = M0Config(quadrature_nodes=11, optimizer_maxiter=400)
    hierarchical = HierarchicalStudentTM0(config, include_study_variance=True).fit(frame)
    h = hierarchical.predict_distribution(frame.iloc[:5])
    assert np.all(h["between_study_sd"] > 0)
    assert np.all(h["upper_80"] > h["prediction"])
    assert np.all(h["lower_80"] < h["prediction"])
    # The hierarchical fit may assign less residual noise than a non-hierarchical
    # fit, so only the within-model variance decomposition is a valid invariant.
    residual_variance = h["residual_scale"] ** 2 * config.degrees_of_freedom / (config.degrees_of_freedom - 2)
    assert np.all(residual_variance + h["between_study_sd"] ** 2 > residual_variance)


def test_fit_and_prediction_are_reproducible() -> None:
    frame = synthetic_frame()
    config = M0Config(quadrature_nodes=9, optimizer_maxiter=350)
    one = HierarchicalStudentTM0(config).fit(frame)
    two = HierarchicalStudentTM0(config).fit(frame)
    np.testing.assert_allclose(one.predict_location(frame), two.predict_location(frame), rtol=0, atol=1e-10)
    assert json.dumps(one.parameters(), sort_keys=True) == json.dumps(two.parameters(), sort_keys=True)


def test_leakage_contract_rejects_study_compound_and_temporal_overlap() -> None:
    frame = synthetic_frame().iloc[:8].copy().reset_index(drop=True)
    with pytest.raises(AssertionError, match="Study leakage"):
        _assert_fold_safe(frame, np.asarray([0, 1]), np.asarray([2]), "study_grouped")
    frame.loc[3, "molecule_id"] = frame.loc[0, "molecule_id"]
    with pytest.raises(AssertionError, match="Compound leakage"):
        _assert_fold_safe(frame, np.asarray([0]), np.asarray([3]), "study_grouped")
    temporal = synthetic_frame().reset_index(drop=True)
    with pytest.raises(AssertionError, match="Future-information leakage"):
        _assert_fold_safe(temporal, np.asarray([30]), np.asarray([0]), "forward_temporal")


def test_metric_calculation_preserves_natural_product_scope() -> None:
    frame = pd.DataFrame({
        "validation_scheme": ["study_grouped"] * 4,
        "method": ["phytotrust_m0"] * 4,
        "heldout_study_key": ["A", "B", "C", "D"],
        "heldout_median_pactivity": [5.0, 6.0, 7.0, 8.0],
        "prediction": [5.1, 5.9, 6.8, 8.2],
        "nll": [0.2, 0.2, 0.3, 0.3],
        "lower_80": [4.0, 5.0, 6.0, 7.0], "upper_80": [6.0, 7.0, 8.0, 9.0],
        "lower_95": [3.0, 4.0, 5.0, 6.0], "upper_95": [7.0, 8.0, 9.0, 10.0],
        "is_verified_natural_product": [True, False, True, False],
    })
    _, metrics = evaluate_predictions(frame)
    assert metrics["study_grouped"]["all"]["phytotrust_m0"]["n"] == 4
    assert metrics["study_grouped"]["natural_products"]["phytotrust_m0"]["n"] == 2


def test_study_cluster_bootstrap_is_reproducible() -> None:
    rows = []
    methods = ("global_mean", "median_potency", "nonhierarchical_student_t", "phytotrust_m0")
    for scheme in ("study_grouped", "forward_temporal"):
        for index in range(8):
            observed = 5.0 + index * 0.1
            for method_index, method in enumerate(methods):
                prediction = observed + (method_index - 2) * 0.04
                rows.append({
                    "validation_scheme": scheme, "source_index": index,
                    "heldout_study_key": f"S{index // 2}", "method": method,
                    "heldout_median_pactivity": observed, "prediction": prediction,
                    "nll": 0.3 + abs(prediction - observed),
                    "lower_80": prediction - 0.5, "upper_80": prediction + 0.5,
                    "is_verified_natural_product": index % 2 == 0,
                })
    predictions = pd.DataFrame(rows)
    first, first_nested = cluster_bootstrap_differences(predictions, iterations=25, seed=91)
    second, second_nested = cluster_bootstrap_differences(predictions, iterations=25, seed=91)
    pd.testing.assert_frame_equal(first, second)
    assert first_nested == second_nested
    assert "point" in first_nested["study_grouped"]["all"]["median_potency"]
    assert "ci95" in first_nested["study_grouped"]["all"]["median_potency"]
