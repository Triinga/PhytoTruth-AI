from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from app.phytotrust_m0 import M0Config
from app.phytotrust_m0_connected import (
    aggregate_compound_study,
    assert_fold_provenance,
    construct_test_events,
    construct_training_events,
    cross_validate,
    evaluate_predictions,
    load_connected_development,
)


def measurement_frame() -> pd.DataFrame:
    rows = []
    activity = 0
    for study_index, study in enumerate(("D1", "D2", "D3", "D4", "D5")):
        for compound_index in range(6):
            activity += 1
            rows.append({
                "activity_id": str(activity), "target_chembl_id": "CHEMBL222",
                "molecule_id": f"M{compound_index}", "assay_id": f"A{study_index}",
                "document_id": study, "document_year": 2008 + study_index,
                "pchembl_value": 5.0 + compound_index * 0.2 + study_index * 0.1,
                "partition": "development", "assay_context_key": "CTX",
                "component_id": "C1", "is_verified_natural_product": compound_index == 0,
                "natural_product_sources": "NPASS" if compound_index == 0 else "",
            })
    return pd.DataFrame(rows)


def test_loader_fails_closed_before_reading_quarantined_file(tmp_path: Path) -> None:
    forbidden = tmp_path / "QUARANTINED_future_evaluation_labels.csv"
    forbidden.write_text("not,a,dataset\n", encoding="utf-8")
    with pytest.raises(ValueError, match="evaluation artifacts are forbidden"):
        load_connected_development(forbidden)


def test_loader_accepts_only_development_partition_and_preserves_hash(tmp_path: Path) -> None:
    path = tmp_path / "development_measurements.csv"
    measurement_frame().to_csv(path, index=False)
    frame, manifest = load_connected_development(path)
    assert len(frame) == 30
    assert manifest["quarantined_files_opened"] is False
    assert len(manifest["source_sha256"]) == 64
    compromised = measurement_frame()
    compromised.loc[0, "partition"] = "future_evaluation"
    compromised.to_csv(path, index=False)
    with pytest.raises(AssertionError, match="Non-development"):
        load_connected_development(path)


def test_aggregation_and_training_events_preserve_publication_provenance() -> None:
    edges = aggregate_compound_study(measurement_frame())
    assert len(edges) == 30
    events = construct_training_events(edges)
    assert len(events) == 30
    for row in events.itertuples(index=False):
        assert row.heldout_document_ids not in row.prior_document_ids.split(";")
        assert not set(row.heldout_activity_ids.split(";")) & set(row.prior_activity_ids.split(";"))


def test_test_event_uses_only_allowed_training_measurements() -> None:
    edges = aggregate_compound_study(measurement_frame())
    train = edges[edges.document_id != "D5"]
    test = edges[edges.document_id == "D5"]
    events = construct_test_events(train, test, history_mode="strictly_earlier_publications")
    assert len(events) == 6
    assert set(events.prior_study_count) == {4}
    assert set(events.prior_document_ids.str.split(";").explode()) == {"D1", "D2", "D3", "D4"}
    assert_fold_provenance(train, events, temporal=True)


def test_temporal_leakage_check_rejects_later_training_publication() -> None:
    edges = aggregate_compound_study(measurement_frame())
    train = edges[edges.document_id != "D3"]
    test = edges[edges.document_id == "D3"]
    events = construct_test_events(train, test, history_mode="invalid_future_history")
    with pytest.raises(AssertionError, match="Forward-temporal"):
        assert_fold_provenance(train, events, temporal=True)


def test_connected_cross_validation_is_reproducible_and_emits_80_90_intervals() -> None:
    edges = aggregate_compound_study(measurement_frame())
    config = M0Config(quadrature_nodes=7, optimizer_maxiter=250)
    first, first_parameters, _ = cross_validate(edges, config)
    second, second_parameters, _ = cross_validate(edges, config)
    sort_columns = ["validation_scheme", "method", "prediction_id"]
    pd.testing.assert_frame_equal(
        first.sort_values(sort_columns).reset_index(drop=True),
        second.sort_values(sort_columns).reset_index(drop=True),
    )
    pd.testing.assert_frame_equal(first_parameters, second_parameters)
    assert {"lower_80", "upper_80", "lower_90", "upper_90"} <= set(first.columns)
    assert not first[first.method == "phytotrust_m0_connected"].empty
    assert np.isfinite(first.nll).all()
    metric_table, nested = evaluate_predictions(first)
    assert not metric_table.empty
    assert nested["study_grouped"]["all"]["phytotrust_m0_connected"]["n"] == 30
