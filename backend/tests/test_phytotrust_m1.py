from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from app.assay_protocol import extract_protocol_record
from app.phytotrust_m1 import (
    ALL_TABULAR_FEATURES,
    EvidenceConditionedM1,
    M1Config,
    construct_feature_events,
    enrich_edges,
    load_protocol_features,
)


def measurements() -> pd.DataFrame:
    rows = []
    activity = 0
    for study_index, study in enumerate(("D1", "D2", "D3", "D4", "D5")):
        for compound_index in range(6):
            activity += 1
            rows.append({
                "activity_id": str(activity), "target_chembl_id": "CHEMBL222",
                "target_name": "NET", "target_type": "SINGLE PROTEIN",
                "target_organism": "Homo sapiens", "molecule_id": f"M{compound_index}",
                "standardized_smiles": "CCO" if compound_index % 2 else "c1ccccc1O",
                "assay_id": f"A{study_index}", "document_id": study,
                "document_year": 2008 + study_index,
                "pchembl_value": 5.0 + compound_index * 0.2 + study_index * 0.12,
                "assay_type": "B", "bao_format": "BAO_0000219",
                "bao_label": "cell-based format",
                "assay_description": "Radioligand binding in HEK293 cells",
                "assay_context_key": "CTX", "component_id": "C1",
                "partition": "development",
                "is_verified_natural_product": compound_index == 0,
                "natural_product_sources": "NPASS" if compound_index == 0 else "",
            })
    return pd.DataFrame(rows)


def protocol_frame() -> pd.DataFrame:
    rows = []
    for study_index in range(5):
        rows.append(extract_protocol_record({"record": {
            "assay_chembl_id": f"A{study_index}", "target_chembl_id": "CHEMBL222",
            "document_chembl_id": f"D{study_index + 1}", "assay_type": "B",
            "assay_organism": "Homo sapiens", "bao_format": "BAO_0000219",
            "bao_label": "cell-based format",
            "description": "Radioligand binding in HEK293 cells at 25 C for 30 minutes",
        } }))
    return pd.DataFrame(rows)


def test_protocol_loader_fails_closed_for_quarantined_path(tmp_path: Path) -> None:
    path = tmp_path / "QUARANTINED_protocol.csv"
    path.write_text("assay_chembl_id\nA0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="forbidden"):
        load_protocol_features(path, {"A0"})


def test_protocol_join_preserves_real_assay_provenance_and_missingness() -> None:
    edges = enrich_edges(measurements(), protocol_frame())
    assert len(edges) == 30
    assert set(edges.protocol_source_assay_ids) == {"A0", "A1", "A2", "A3", "A4"}
    assert edges.protocol_missing_fraction.between(0, 1).all()
    assert (edges.protocol_category_count > 0).all()


def test_heldout_protocol_never_changes_historical_features() -> None:
    edges = enrich_edges(measurements(), protocol_frame())
    train, test = edges[edges.document_id != "D5"], edges[edges.document_id == "D5"].copy()
    first = construct_feature_events(train, test, history_mode="strictly_earlier_publications")
    test["protocol_missing_fraction"] = 0.0
    test["protocol_category_count"] = 999
    test["protocol_tokens"] = "future:leak"
    second = construct_feature_events(train, test, history_mode="strictly_earlier_publications")
    pd.testing.assert_frame_equal(first, second)


def synthetic_events(seed: int = 9) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for study_index in range(10):
        study_shift = rng.normal(0, 0.3)
        for compound_index in range(4):
            features = {
                "history_mean_minus_median": rng.normal(0, 0.05),
                "history_mad": rng.uniform(0, 0.6),
                "log1p_prior_studies": rng.uniform(0.69, 1.4),
                "molecular_weight_scaled": rng.uniform(0.3, 0.9),
                "logp_scaled": rng.uniform(0.1, 1.1),
                "inverse_sqrt_prior_studies": rng.uniform(0.5, 1.0),
                "protocol_missing_fraction": rng.uniform(0.4, 1.0),
                "protocol_category_diversity": rng.uniform(0, 3),
                "history_year_span_scaled": rng.uniform(0, 1),
            }
            prior = 6.0 + rng.normal(0, 0.4)
            scale = 0.1 + 0.35 * features["protocol_missing_fraction"]
            rows.append({
                "target_chembl_id": "CHEMBL222", "molecule_id": f"M{compound_index}",
                "heldout_study_key": f"S{study_index}", "heldout_year": 2000 + study_index,
                "prior_median_pactivity": prior,
                "heldout_median_pactivity": prior + study_shift + rng.standard_t(4) * scale,
                **features,
            })
    return pd.DataFrame(rows)


def test_m1_is_reproducible_and_returns_uncertainty_decomposition() -> None:
    frame = synthetic_events()
    config = M1Config(quadrature_nodes=9, optimizer_maxiter=400)
    first = EvidenceConditionedM1(config).fit(frame)
    second = EvidenceConditionedM1(config).fit(frame)
    assert first.parameters() == second.parameters()
    distribution = first.predict_distribution(frame.iloc[:6])
    assert {"residual_variance", "between_study_variance", "lower_80", "upper_90"} <= set(distribution)
    assert np.all(distribution["residual_variance"] > 0)
    assert np.all(distribution["between_study_variance"] > 0)


def test_no_study_ablation_has_zero_between_study_uncertainty() -> None:
    frame = synthetic_events()
    model = EvidenceConditionedM1(
        M1Config(quadrature_nodes=7, optimizer_maxiter=350), include_study_variance=False,
    ).fit(frame)
    distribution = model.predict_distribution(frame.iloc[:5])
    assert np.all(distribution["between_study_sd"] == 0)
    assert np.all(distribution["between_study_variance"] == 0)


def test_m1_rejects_missing_prediction_time_features() -> None:
    frame = synthetic_events().drop(columns=[ALL_TABULAR_FEATURES[-1]])
    with pytest.raises(ValueError, match="Missing M1 columns"):
        EvidenceConditionedM1().fit(frame)
