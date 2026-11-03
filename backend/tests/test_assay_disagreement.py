import pandas as pd

from app.assay_disagreement import (
    ALL_NUMERIC_FEATURES,
    MODEL_FEATURES,
    _case_from_records,
    clustered_metric_intervals,
    clustered_model_comparison,
    fit_models,
)


def activity(activity_id, year, value, document, assay):
    return {
        "activity_id": activity_id,
        "molecule_chembl_id": "CHEMBL_M1",
        "canonical_smiles": "CCO",
        "computed_pactivity": value,
        "document_chembl_id": document,
        "assay_chembl_id": assay,
        "raw_source": {
            "document_year": year,
            "document_journal": "Journal",
            "assay_description": "exact IC50 assay",
        },
    }


def test_case_uses_independent_future_study_and_defines_tenfold_disagreement():
    records = [
        activity(1, 2017, 6.0, "DOC1", "ASSAY1"),
        activity(2, 2018, 6.2, "DOC2", "ASSAY2"),
        activity(3, 2021, 7.2, "DOC3", "ASSAY3"),
    ]
    case, reason = _case_from_records(
        {"target_chembl_id": "CHEMBL1", "target_name": "Target"},
        "CHEMBL_M1", records, {"CHEMBL_M1": [{"source": "NPASS"}]},
    )
    assert reason is None
    assert case["disagreement_label"] == 1
    assert case["partition"] == "evaluation"
    assert case["prior_document_ids"] == "DOC1;DOC2"
    assert case["heldout_document_ids"] == "DOC3"
    assert case["is_verified_natural_product"] is True


def test_case_rejects_assay_reuse_across_time():
    records = [
        activity(1, 2018, 6.0, "DOC1", "ASSAY1"),
        activity(2, 2021, 7.2, "DOC2", "ASSAY1"),
    ]
    case, reason = _case_from_records(
        {"target_chembl_id": "CHEMBL1", "target_name": "Target"},
        "CHEMBL_M1", records, {},
    )
    assert case is None
    assert reason == "assay_overlap"


def test_feature_contract_cannot_read_heldout_outcome_or_metadata():
    all_features = set(ALL_NUMERIC_FEATURES) | {"target_chembl_id"}
    assert not any("heldout" in name for name in all_features)
    assert "absolute_pactivity_shift" not in all_features
    assert "disagreement_label" not in all_features
    assert "is_verified_natural_product" not in all_features
    assert set(MODEL_FEATURES["full"]) == set(ALL_NUMERIC_FEATURES)


def test_clustered_inference_is_reproducible_and_study_aware():
    frame = pd.DataFrame({
        "target_chembl_id": ["T1"] * 6 + ["T2"] * 6,
        "heldout_study_key": ["A"] * 3 + ["B"] * 3 + ["C"] * 3 + ["D"] * 3,
        "disagreement_label": [0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
        "candidate": [0.1, 0.2, 0.1, 0.8, 0.7, 0.9, 0.2, 0.3, 0.8, 0.7, 0.6, 0.4],
        "reference": [0.5] * 12,
    })
    first = clustered_metric_intervals(frame, "candidate", iterations=50, seed=7)
    second = clustered_metric_intervals(frame, "candidate", iterations=50, seed=7)
    assert first == second
    comparison = clustered_model_comparison(
        frame, "candidate", "reference", iterations=50, seed=7
    )
    assert comparison["auroc"]["difference"] > 0


def test_model_predictions_are_deterministic():
    rows = []
    for index in range(40):
        row = {name: float((index + offset) % 7) for offset, name in enumerate(ALL_NUMERIC_FEATURES)}
        row.update({
            "target_chembl_id": f"T{index % 2}",
            "disagreement_label": int(index % 5 == 0),
            "molecule_id": f"M{index}",
        })
        rows.append(row)
    development = pd.DataFrame(rows[:30])
    evaluation = pd.DataFrame(rows[30:])
    _, first = fit_models(development, evaluation)
    _, second = fit_models(development, evaluation)
    assert first.filter(like="probability_").equals(second.filter(like="probability_"))
