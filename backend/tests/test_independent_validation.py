import json

import numpy as np
import pandas as pd
import pytest

from app.independent_validation import (
    assert_independence,
    event_overlap_reasons,
    evaluate_frozen_predictions,
    score_frozen_models,
    strict_independent_filter,
    validation_measurements,
)


def _event(**updates):
    value = {
        "target_chembl_id": "NEW-TARGET",
        "target_name": "New target",
        "molecule_id": "NEW-MOL",
        "heldout_year": 2024,
        "partition": "evaluation",
        "heldout_study_key": "NEW-DOC",
        "prior_activity_ids": "NEW-A1",
        "heldout_activity_ids": "NEW-A2",
        "prior_assay_ids": "NEW-S1",
        "heldout_assay_ids": "NEW-S2",
        "prior_document_ids": "NEW-D1",
        "heldout_document_ids": "NEW-D2",
        "disagreement_label": 0,
        "is_verified_natural_product": False,
        "feature": 1.0,
    }
    value.update(updates)
    return value


def _ledger():
    return {
        "targets": {"OLD-TARGET"},
        "molecules": {"OLD-MOL"},
        "activities": {"OLD-A"},
        "assays": {"OLD-S"},
        "documents": {"OLD-DOC"},
    }


@pytest.mark.parametrize(
    "updates,reason",
    [
        ({"target_chembl_id": "OLD-TARGET"}, "previous_target"),
        ({"molecule_id": "OLD-MOL"}, "previous_compound"),
        ({"prior_activity_ids": "OLD-A"}, "previous_activity"),
        ({"heldout_assay_ids": "OLD-S"}, "previous_assay"),
        ({"heldout_document_ids": "OLD-DOC"}, "previous_publication"),
    ],
)
def test_every_prior_exposure_type_is_detected(updates, reason):
    assert reason in event_overlap_reasons(_event(**updates), _ledger())


def test_strict_filter_enforces_temporal_and_within_event_leakage():
    candidates = pd.DataFrame([
        _event(molecule_id="GOOD"),
        _event(molecule_id="OLD-MOL"),
        _event(molecule_id="EARLY", heldout_year=2019, partition="development"),
        _event(molecule_id="ASSAY", heldout_assay_ids="NEW-S1"),
        _event(molecule_id="DOC", heldout_document_ids="NEW-D1"),
    ])
    accepted, rejected, counts = strict_independent_filter(candidates, _ledger(), minimum_heldout_year=2020)
    assert accepted.molecule_id.tolist() == ["GOOD"]
    assert counts["previous_compound"] == 1
    assert counts["heldout_before_frozen_temporal_boundary"] == 1
    assert counts["within_event_assay_overlap"] == 1
    assert counts["within_event_publication_overlap"] == 1
    assert len(rejected) == 4
    assert_independence(accepted, _ledger(), minimum_heldout_year=2020)


def test_compounds_repeated_across_new_targets_are_all_excluded():
    candidates = pd.DataFrame([
        _event(target_chembl_id="T1", molecule_id="DUP"),
        _event(target_chembl_id="T2", molecule_id="DUP"),
        _event(target_chembl_id="T3", molecule_id="UNIQUE"),
    ])
    accepted, rejected, counts = strict_independent_filter(candidates, _ledger(), minimum_heldout_year=2020)
    assert accepted.molecule_id.tolist() == ["UNIQUE"]
    assert counts["compound_repeated_across_new_targets"] == 2
    assert len(rejected) == 2


class FrozenDummyModel:
    def predict_proba(self, frame):
        probabilities = np.asarray(frame["feature"], dtype=float)
        return np.column_stack([1 - probabilities, probabilities])


def test_frozen_scoring_uses_artifact_without_fitting():
    frame = pd.DataFrame([_event(molecule_id="M1", feature=0.2), _event(molecule_id="M2", feature=0.8)])
    full = {
        "model": FrozenDummyModel(), "features": ["feature"],
        "validation_status": "not_supported_for_reliable_prediction",
    }
    selected = {
        "constant_prevalence": 0.25,
        "validation_status": "not_supported_for_reliable_prediction",
    }
    scored = score_frozen_models(frame, full, selected)
    assert scored.probability_frozen_full.tolist() == pytest.approx([0.2, 0.8])
    assert scored.probability_frozen_selected.tolist() == pytest.approx([0.25, 0.25])
    assert not hasattr(full["model"], "fit_called")


def test_clustered_evaluation_is_seed_reproducible():
    rows = []
    for target in ("T1", "T2"):
        for study in range(3):
            for index in range(2):
                label = (study + index) % 2
                rows.append({
                    **_event(target_chembl_id=target, molecule_id=f"{target}-{study}-{index}"),
                    "heldout_study_key": f"{target}-S{study}",
                    "disagreement_label": label,
                    "probability_frozen_full": 0.8 if label else 0.2,
                    "probability_frozen_selected": 0.25,
                })
    frame = pd.DataFrame(rows)
    first = evaluate_frozen_predictions(frame, iterations=30, seed=17)
    second = evaluate_frozen_predictions(frame, iterations=30, seed=17)
    assert first == second
    assert first["overall"]["frozen_full"]["point"]["auroc"] == pytest.approx(1.0)


def test_used_measurements_preserve_event_and_role_provenance(tmp_path):
    activities = tmp_path / "activities.jsonl"
    records = [
        {"activity_id": "A1", "molecule_chembl_id": "M1", "assay_chembl_id": "S1",
         "document_chembl_id": "D1", "computed_pactivity": 6.0, "raw_source": {"document_year": 2020}},
        {"activity_id": "A2", "molecule_chembl_id": "M1", "assay_chembl_id": "S2",
         "document_chembl_id": "D2", "computed_pactivity": 7.2, "raw_source": {"document_year": 2024}},
    ]
    activities.write_text("".join(json.dumps(row) + "\n" for row in records), encoding="utf-8")
    natural = tmp_path / "natural.jsonl"
    natural.write_text("", encoding="utf-8")
    specs = tmp_path / "specs.json"
    specs.write_text(json.dumps([{
        "target_chembl_id": "T1", "target_name": "Target", "organism": "Homo sapiens",
        "target_type": "SINGLE PROTEIN", "activities": str(activities), "natural_products": str(natural),
    }]), encoding="utf-8")
    events = pd.DataFrame([{**_event(), "validation_event_id": "E1", "prior_activity_ids": "A1", "heldout_activity_ids": "A2"}])
    result = validation_measurements(events, specs)
    assert result.activity_id.tolist() == ["A2", "A1"] or set(result.activity_id) == {"A1", "A2"}
    assert dict(zip(result.activity_id, result.evidence_role)) == {"A1": "prior", "A2": "heldout"}
    assert set(result.validation_event_id) == {"E1"}
