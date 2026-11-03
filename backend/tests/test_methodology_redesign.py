import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from app.methodology_redesign import (
    FULL_FEATURES,
    _fit_fold,
    assert_feature_contract,
    load_development_dataset,
    make_splits,
)


ROOT = Path(__file__).resolve().parents[2]
DATASET = ROOT / "experiments" / "assay_protocol_enrichment" / "enriched_assay_disagreement_events.csv"


@pytest.fixture(scope="module")
def development() -> pd.DataFrame:
    frame, manifest = load_development_dataset(DATASET)
    assert manifest["excluded_locked_evaluation_rows"] > 0
    return frame


def test_loader_excludes_locked_evaluation_and_preserves_threshold(development: pd.DataFrame) -> None:
    assert set(development.partition) == {"development"}
    assert development.heldout_year.max() <= 2019
    assert np.array_equal(
        (development.absolute_pactivity_shift >= 1.0).astype(int),
        development.disagreement_label,
    )


def test_feature_contract_is_prior_only() -> None:
    assert_feature_contract(FULL_FEATURES)
    with pytest.raises(ValueError, match="forbidden"):
        assert_feature_contract(["heldout_parameter_reported_fraction"])
    with pytest.raises(ValueError, match="forbidden"):
        assert_feature_contract(["protocol_organism_comparison"])


def test_all_splits_prevent_compound_and_required_group_leakage(development: pd.DataFrame) -> None:
    splits = make_splits(development)
    assert set(splits) == {"study_grouped", "forward_temporal", "leave_target_out"}
    for scheme, folds in splits.items():
        for train, test in folds:
            assert set(development.iloc[train].molecule_id).isdisjoint(development.iloc[test].molecule_id)
            if scheme == "study_grouped":
                assert set(development.iloc[train].heldout_study_key).isdisjoint(development.iloc[test].heldout_study_key)
            if scheme == "forward_temporal":
                assert development.iloc[train].heldout_year.max() < development.iloc[test].heldout_year.min()
            if scheme == "leave_target_out":
                assert set(development.iloc[train].target_chembl_id).isdisjoint(development.iloc[test].target_chembl_id)


def test_fold_predictions_are_deterministic_and_intervals_ordered(development: pd.DataFrame) -> None:
    train, test = make_splits(development)["study_grouped"][0]
    first = _fit_fold(development.iloc[train], development.iloc[test], seed=2026)
    second = _fit_fold(development.iloc[train], development.iloc[test], seed=2026)
    assert set(first) == set(second)
    for name in first:
        np.testing.assert_allclose(first[name], second[name], rtol=0, atol=0)
        assert np.isfinite(first[name]).all()
    for model in ("rf_history", "rf_protocol"):
        assert np.all(first[f"lower_{model}"] <= first[f"magnitude_{model}"])
        assert np.all(first[f"magnitude_{model}"] <= first[f"upper_{model}"])
        assert np.all((first[f"risk_{model}"] >= 0) & (first[f"risk_{model}"] <= 1))


def test_executed_results_record_locked_data_exclusion() -> None:
    results_path = ROOT / "experiments" / "assay_disagreement_methodology_redesign" / "results.json"
    results = json.loads(results_path.read_text(encoding="utf-8"))
    assert results["independent_evaluation_datasets_loaded"] is False
    assert results["dataset"]["development_rows"] == 212
    assert results["dataset"]["excluded_locked_evaluation_rows"] == 319
    assert results["readiness"]["status"] == "not_ready_for_new_independent_test"
