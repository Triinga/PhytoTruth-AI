import numpy as np
import pandas as pd

from app.assay_disagreement_generalization import (
    COMPACT_FEATURES,
    cross_fitted_calibration,
    leave_one_study_out_indices,
    study_balanced_metrics,
)


def test_leave_one_study_out_indices_are_group_disjoint_and_complete():
    studies = np.asarray(["A", "A", "B", "C", "C", "D"])
    splits = leave_one_study_out_indices(studies)
    seen = []
    assert len(splits) == 4
    for train, test in splits:
        assert set(studies[train]).isdisjoint(set(studies[test]))
        assert len(set(studies[test])) == 1
        seen.extend(test.tolist())
    assert sorted(seen) == list(range(len(studies)))


def test_cross_fitted_calibration_is_deterministic_and_finite():
    raw = np.asarray([0.1, 0.2, 0.8, 0.7, 0.3, 0.6, 0.15, 0.85, 0.4, 0.65])
    labels = np.asarray([0, 0, 1, 1, 0, 1, 0, 1, 0, 1])
    studies = np.asarray(["A", "A", "B", "B", "C", "C", "D", "D", "E", "E"])
    first, _ = cross_fitted_calibration(raw, labels, studies, seed=11)
    second, _ = cross_fitted_calibration(raw, labels, studies, seed=11)
    assert np.allclose(first, second)
    assert np.all((first > 0) & (first < 1))


def test_study_balanced_metric_prevents_large_study_domination():
    labels = [1] * 10 + [0]
    probabilities = [0.9] * 10 + [0.9]
    studies = ["large"] * 10 + ["small"]
    result = study_balanced_metrics(labels, probabilities, studies)
    expected = (0.01 + 0.81) / 2
    assert abs(result["brier"] - expected) < 1e-12


def test_enhanced_feature_contract_has_no_future_or_np_fields():
    prohibited_tokens = ("heldout", "future", "disagreement_label", "natural_product")
    assert not any(token in feature for feature in COMPACT_FEATURES for token in prohibited_tokens)
    assert "prior_recent_shift" in COMPACT_FEATURES
    assert "aromatic_atom_fraction" in COMPACT_FEATURES
