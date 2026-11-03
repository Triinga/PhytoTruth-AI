import json

import pandas as pd
import pytest

from app.assay_protocol import (
    association_analysis,
    enrich_events,
    extract_protocol_record,
    load_and_extract,
    missingness_association_analysis,
    run_enrichment,
)


def _envelope(assay_id, description=None, **fields):
    record = {
        "assay_chembl_id": assay_id,
        "target_chembl_id": "CHEMBL1",
        "document_chembl_id": f"DOC-{assay_id}",
        "assay_type": "B",
        "assay_type_description": "Binding",
        "description": description,
        "assay_parameters": [],
        "bao_format": "BAO_0000357",
        "bao_label": "single protein format",
        **fields,
    }
    return {
        "source": "ChEMBL REST API",
        "source_endpoint": f"https://example.test/{assay_id}",
        "retrieved_at": "2026-01-01T00:00:00+00:00",
        "record": record,
    }


def test_literal_protocol_extraction_preserves_evidence_and_missingness():
    row = extract_protocol_record(_envelope(
        "CHEMBL-A",
        "Radioligand binding using recombinant membranes; substrate ATP at pH 7.4, incubated 30 min at 25 C.",
    ))
    assert row["detection_categories"] == "radiometric"
    assert "recombinant" in row["biological_system_text_categories"]
    assert row["ph_mentions"] == "pH 7.4"
    assert "30 min" in row["duration_mentions"]
    assert row["temperature_mentions"] == "25 C"
    evidence = json.loads(row["extraction_evidence_json"])
    assert evidence["source_field"] == "description"
    assert "substrate ATP" in row["substrate_mentions_json"]

    missing = extract_protocol_record(_envelope("CHEMBL-B"))
    assert missing["description_status"] == "missing"
    assert missing["detection_categories"] == ""
    assert missing["text_extraction_status"] == "no_source_description"


def test_rule_engine_does_not_infer_unmentioned_detection_or_conditions():
    row = extract_protocol_record(_envelope("CHEMBL-A", "Inhibition of kinase activity"))
    assert row["methodology_categories"] == "enzyme_activity"
    assert row["detection_categories"] == ""
    assert row["temperature_mentions"] == ""
    assert row["ph_mentions"] == ""
    assert row["substrate_mentions_json"] == "[]"


def test_duplicate_assay_metadata_is_rejected(tmp_path):
    path = tmp_path / "raw.jsonl"
    envelope = _envelope("CHEMBL-A", "Binding")
    path.write_text(json.dumps(envelope) + "\n" + json.dumps(envelope) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate assay metadata"):
        load_and_extract(path)


def test_event_enrichment_labels_heldout_protocol_as_post_outcome():
    assays = pd.DataFrame([
        extract_protocol_record(_envelope("CHEMBL-A", "Radioligand binding", assay_organism="Homo sapiens")),
        extract_protocol_record(_envelope("CHEMBL-B", "Fluorescence binding", assay_organism="Rattus norvegicus")),
    ])
    events = pd.DataFrame([{
        "target_chembl_id": "CHEMBL1", "molecule_id": "M1", "heldout_study_key": "DOC2",
        "prior_assay_ids": "CHEMBL-A", "heldout_assay_ids": "CHEMBL-B", "partition": "evaluation",
        "disagreement_label": 1, "is_verified_natural_product": False,
    }])
    result = enrich_events(events, assays).iloc[0]
    assert result.protocol_organism_changed == 1
    assert result.protocol_detection_changed == 1
    assert result.protocol_prior_feature_timing == "pre_prediction_available"
    assert result.protocol_heldout_feature_timing == "post_outcome_explanatory_only"
    assert result.protocol_comparison_feature_timing == "post_outcome_explanatory_only"
    assert result.prior_detection_reported == 1
    assert result.heldout_detection_reported == 1


def test_unknown_protocol_comparison_remains_unknown():
    assays = pd.DataFrame([
        extract_protocol_record(_envelope("CHEMBL-A", None)),
        extract_protocol_record(_envelope("CHEMBL-B", None)),
    ])
    events = pd.DataFrame([{
        "target_chembl_id": "CHEMBL1", "molecule_id": "M1", "heldout_study_key": "DOC2",
        "prior_assay_ids": "CHEMBL-A", "heldout_assay_ids": "CHEMBL-B", "partition": "evaluation",
        "disagreement_label": 0, "is_verified_natural_product": False,
    }])
    result = enrich_events(events, assays).iloc[0]
    assert result.protocol_detection_comparison == "unknown"
    assert pd.isna(result.protocol_detection_changed)


def test_association_analysis_is_seed_reproducible():
    rows = []
    for study, changed, outcomes in [("A", 0, [0, 0]), ("B", 1, [1, 0]), ("C", 1, [1, 1])]:
        for outcome in outcomes:
            rows.append({
                "partition": "evaluation", "heldout_study_key": study,
                "disagreement_label": outcome, "protocol_detection_changed": changed,
                "is_verified_natural_product": False,
            })
    frame = pd.DataFrame(rows)
    first = association_analysis(frame, iterations=30, seed=7)
    second = association_analysis(frame, iterations=30, seed=7)
    assert first == second
    assert first[0]["feature_timing"] == "post_outcome_explanatory_only"


def test_prior_missingness_analysis_is_prediction_time_only():
    frame = pd.DataFrame([
        {"partition": "evaluation", "heldout_study_key": f"S{i // 2}", "disagreement_label": i % 2,
         "prior_detection_reported": int(i >= 2), "is_verified_natural_product": False}
        for i in range(6)
    ])
    result = missingness_association_analysis(frame, iterations=20, seed=9)
    assert result[0]["feature_timing"] == "pre_prediction_available"
    assert "risk_difference_reported_minus_missing" in result[0]


def test_full_enrichment_outputs_are_reproducible(tmp_path):
    events = pd.DataFrame([
        {"target_chembl_id": "CHEMBL1", "molecule_id": f"M{i}", "heldout_study_key": f"DOC{i % 3}",
         "prior_assay_ids": "CHEMBL-A", "heldout_assay_ids": "CHEMBL-B", "partition": "evaluation",
         "disagreement_label": i % 2, "is_verified_natural_product": i == 0}
        for i in range(6)
    ])
    event_path = tmp_path / "events.csv"
    events.to_csv(event_path, index=False)
    raw_path = tmp_path / "raw.jsonl"
    raw_path.write_text(
        json.dumps(_envelope("CHEMBL-A", "Radioligand binding")) + "\n" +
        json.dumps(_envelope("CHEMBL-B", "Fluorescence binding")) + "\n",
        encoding="utf-8",
    )
    first = run_enrichment(event_path, raw_path, tmp_path / "one", bootstrap_iterations=20, seed=11)
    second = run_enrichment(event_path, raw_path, tmp_path / "two", bootstrap_iterations=20, seed=11)
    assert first["coverage"] == second["coverage"]
    assert (tmp_path / "one" / "enriched_assay_disagreement_events.csv").read_bytes() == (
        tmp_path / "two" / "enriched_assay_disagreement_events.csv"
    ).read_bytes()
