import json
from pathlib import Path

import pandas as pd

from app.dataset_expansion import (
    _event_identifier,
    classify_protocol_explanation,
    merge_ledgers,
)


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "experiments" / "assay_disagreement_expansion"


def _ids(value):
    return {item for item in str(value).split(";") if item}


def test_exposure_ledgers_are_unioned_without_loss():
    first = {key: {f"a-{key}"} for key in ("targets", "molecules", "activities", "assays", "documents")}
    second = {key: {f"b-{key}"} for key in first}
    merged = merge_ledgers(first, second)
    for key in first:
        assert merged[key] == first[key] | second[key]


def test_event_identifier_is_deterministic_and_outcome_blind():
    base = {
        "target_chembl_id": "T1", "molecule_id": "M1",
        "heldout_study_key": "D2", "heldout_year": 2024,
    }
    assert _event_identifier(base) == _event_identifier({**base, "disagreement_label": 1})
    assert _event_identifier(base).startswith("ADE-")


def test_protocol_interpretation_does_not_infer_missing_facts():
    assert classify_protocol_explanation({
        "protocol_methodology_comparison": "different",
        "prior_protocol_availability": "complete",
        "heldout_protocol_availability": "complete",
    }) == "documented_protocol_difference_possible_explanation"
    assert classify_protocol_explanation({
        "protocol_methodology_comparison": "unknown",
        "protocol_detection_comparison": "unknown",
        "prior_protocol_availability": "complete",
        "heldout_protocol_availability": "complete",
    }) == "unresolved_due_to_incomplete_protocol_metadata"
    assert classify_protocol_explanation({
        **{f"protocol_{name}_comparison": "same" for name in ("assay_type", "organism", "bao_format", "methodology")},
        "prior_protocol_availability": "complete",
        "heldout_protocol_availability": "complete",
    }) == "unexplained_with_partial_documented_similarity"


def test_expansion_artifacts_enforce_measurement_and_leakage_invariants():
    events = pd.read_csv(ARTIFACTS / "expansion_events.csv")
    measurements = pd.read_csv(ARTIFACTS / "used_measurements.csv")
    assert len(events) == 94
    assert events.expansion_event_id.is_unique
    assert events.molecule_id.is_unique
    assert (events.heldout_year >= 2020).all()
    assert (events.partition == "evaluation").all()
    assert (events.disagreement_label == (events.absolute_pactivity_shift >= 1.0).astype(int)).all()
    for row in events.itertuples(index=False):
        assert not (_ids(row.prior_assay_ids) & _ids(row.heldout_assay_ids))
        assert not (_ids(row.prior_document_ids) & _ids(row.heldout_document_ids))
    assert set(measurements.expansion_event_id) == set(events.expansion_event_id)
    assert set(measurements.evidence_role) == {"prior", "heldout"}
    assert (measurements.standard_type == "IC50").all()
    assert (measurements.standard_relation == "=").all()
    assert measurements.activity_id.is_unique


def test_full_records_and_natural_product_provenance_are_traceable():
    events = pd.read_csv(ARTIFACTS / "expansion_events.csv")
    full_path = ARTIFACTS / "used_measurements.full.jsonl"
    records = [json.loads(line) for line in full_path.read_text(encoding="utf-8").splitlines() if line]
    assert len(records) == 295
    assert all(item["normalized_record"]["standard_relation"] == "=" for item in records)
    assert all(not item["normalized_record"]["is_censored"] for item in records)
    assert all(not item["normalized_record"].get("data_validity_comment") for item in records)
    provenance = pd.read_csv(ARTIFACTS / "natural_product_provenance.csv")
    natural_ids = set(events.loc[events.is_verified_natural_product, "molecule_id"])
    assert set(provenance.molecule_id) == natural_ids
    assert set(provenance.source) == {"NPASS", "COCONUT"}
    assert provenance.source_id.notna().all()
    assert provenance.standard_inchi_key.notna().all()


def test_manifest_quarantines_labels_from_model_development():
    manifest = json.loads((ARTIFACTS / "final" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["no_model_training_or_evaluation"] is True
    assert manifest["label_quarantine"] == "not available to model development"
