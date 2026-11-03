import pandas as pd
import pytest

from app.expanded_validation import (
    assert_historical_feature_contract,
    audit_exposure,
    clustered_intervals,
    metric_values,
    tie_aware_precision_at_k,
    score_eligible,
)


class _FrozenModel:
    def predict_proba(self, frame):
        import numpy as np
        values = frame["prior_mad"].to_numpy(float)
        return np.column_stack([1 - values, values])


def _event(**updates):
    row = {
        "expansion_event_id": "E1", "target_chembl_id": "T1", "molecule_id": "M1",
        "prior_activity_ids": "A1", "heldout_activity_ids": "A2",
        "prior_assay_ids": "S1", "heldout_assay_ids": "S2",
        "prior_document_ids": "D1", "heldout_document_ids": "D2",
        "heldout_study_key": "D2", "disagreement_label": 1,
        "probability": 0.8,
    }
    row.update(updates)
    return row


def _ledger():
    return {
        "targets": {"OLD-T"}, "molecules": {"OLD-M"}, "activities": {"OLD-A"},
        "assays": {"OLD-S"}, "documents": {"OLD-D"},
    }


def test_tie_aware_precision_is_order_invariant_and_fractional_at_boundary():
    labels = [1, 0, 1, 0]
    scores = [0.8, 0.8, 0.2, 0.1]
    expected = 0.5
    assert tie_aware_precision_at_k(labels, scores, 1) == pytest.approx(expected)
    assert tie_aware_precision_at_k(labels[::-1], scores[::-1], 1) == pytest.approx(expected)
    assert tie_aware_precision_at_k(labels, [0.3] * 4, 2) == pytest.approx(0.5)


def test_future_or_outcome_features_are_rejected():
    assert_historical_feature_contract(["molecular_weight", "prior_mad"])
    with pytest.raises(AssertionError, match="future features"):
        assert_historical_feature_contract(["prior_mad", "heldout_median_pactivity"])
    with pytest.raises(AssertionError, match="future features"):
        assert_historical_feature_contract(["disagreement_label"])


def test_audit_reports_every_exposure_type_and_within_event_overlap():
    frame = pd.DataFrame([
        _event(expansion_event_id="GOOD"),
        _event(expansion_event_id="BAD", target_chembl_id="OLD-T", molecule_id="OLD-M",
               prior_activity_ids="OLD-A", prior_assay_ids="OLD-S", prior_document_ids="OLD-D"),
        _event(expansion_event_id="WITHIN", molecule_id="M3", heldout_assay_ids="S1", heldout_document_ids="D1"),
    ])
    rows, report = audit_exposure(frame, {"prior": _ledger()})
    assert len(rows) == 1
    assert set(rows.iloc[0].reasons.split(";")) == {
        "previous_target", "previous_compound", "previous_activity", "previous_assay", "previous_publication"
    }
    assert report["union"]["within_event_assay_overlap"] == 1
    assert report["union"]["within_event_publication_overlap"] == 1


def test_study_clustered_intervals_are_seed_deterministic():
    rows = []
    for target in ("T1", "T2"):
        for study in range(4):
            for index in range(2):
                label = (study + index) % 2
                rows.append(_event(
                    expansion_event_id=f"{target}-{study}-{index}", target_chembl_id=target,
                    molecule_id=f"M-{target}-{study}-{index}", heldout_study_key=f"{target}-D{study}",
                    disagreement_label=label, probability=0.8 if label else 0.2,
                ))
    frame = pd.DataFrame(rows)
    first = clustered_intervals(frame, "probability", k=5, iterations=40, seed=31)
    second = clustered_intervals(frame, "probability", k=5, iterations=40, seed=31)
    assert first == second
    assert metric_values(frame, "probability", k=5)["auroc"] == pytest.approx(1.0)


def test_constant_score_precision_equals_sample_prevalence():
    frame = pd.DataFrame([
        _event(expansion_event_id=f"E{i}", molecule_id=f"M{i}", disagreement_label=int(i < 3), probability=0.25)
        for i in range(10)
    ])
    assert metric_values(frame, "probability", k=5)["precision_at_5"] == pytest.approx(0.3)


def test_frozen_scoring_is_deterministic_and_does_not_fit():
    frame = pd.DataFrame([
        {**_event(expansion_event_id="E1", molecule_id="M1"), "prior_mad": 0.2,
         "prior_historical_tenfold_fraction": 0.0},
        {**_event(expansion_event_id="E2", molecule_id="M2"), "prior_mad": 0.8,
         "prior_historical_tenfold_fraction": 0.5},
    ])
    full = {"model": _FrozenModel(), "features": ["prior_mad"], "validation_status": "frozen"}
    selected = {"constant_prevalence": 0.25, "validation_status": "frozen"}
    first = score_eligible(frame, full, selected)
    second = score_eligible(frame, full, selected)
    assert first.probability_frozen_full.tolist() == pytest.approx([0.2, 0.8])
    assert first.probability_frozen_full.tolist() == second.probability_frozen_full.tolist()
    assert first.probability_frozen_selected.tolist() == pytest.approx([0.25, 0.25])
