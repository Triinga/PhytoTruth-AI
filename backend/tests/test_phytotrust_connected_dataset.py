from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from app.phytotrust_connected_dataset import (
    assert_partition_integrity,
    connected_components,
    connectivity_metrics,
    load_locked_ledger,
    select_contexts_and_targets,
)


def graph_frame() -> pd.DataFrame:
    rows = []
    edges = {
        "D1": ["M1", "M2", "M3"],
        "D2": ["M1", "M2", "M4"],
        "D3": ["M1", "M3", "M4"],
    }
    activity = 0
    for document, molecules in edges.items():
        for molecule in molecules:
            activity += 1
            rows.append({
                "activity_id": str(activity), "target_chembl_id": "T1",
                "assay_context_key": "CTX", "molecule_id": molecule,
                "document_id": document, "assay_id": f"A{document}",
                "document_year": 2015 + int(document[-1]),
                "pchembl_value": 5.0 + activity / 10,
                "is_verified_natural_product": molecule == "M1",
                "natural_product_sources": "NPASS" if molecule == "M1" else "",
            })
    return pd.DataFrame(rows)


def protocol() -> dict:
    return {
        "target_selection_rule": {
            "minimum_eligible_measurements": 1,
            "minimum_repeated_compounds": 1,
            "minimum_multi_compound_studies": 1,
            "retain": 1,
        },
        "identifiability_gate": {
            "per_target_minimum_multi_event_studies": 1,
            "per_target_minimum_studies_with_at_least_3_compounds": 1,
            "per_target_minimum_repeated_compounds": 1,
            "required_shared_compound_study_cycles": 1,
        },
    }


def test_connectivity_counts_true_bipartite_overlap() -> None:
    metrics = connectivity_metrics(graph_frame())
    assert metrics["graph_edges"] == 9
    assert metrics["repeated_compounds"] == 4
    assert metrics["multi_compound_studies"] == 3
    assert metrics["study_pairs_sharing_at_least_2_compounds"] == 3
    assert metrics["bipartite_four_cycles"] == 3


def test_components_are_deterministic_and_future_component_is_quarantined() -> None:
    frame = graph_frame()
    first_nodes, first_components = connected_components(frame)
    second_nodes, second_components = connected_components(frame.sample(frac=1, random_state=9))
    pd.testing.assert_frame_equal(
        first_nodes.sort_values(["node_type", "node_id"]).reset_index(drop=True),
        second_nodes.sort_values(["node_type", "node_id"]).reset_index(drop=True),
    )
    assert len(first_components) == 1
    assert first_components.iloc[0].partition == "development"
    frame.loc[frame.document_id == "D3", "document_year"] = 2021
    _, future_components = connected_components(frame)
    assert future_components.iloc[0].partition == "future_evaluation"


def test_context_and_target_selection_do_not_depend_on_activity_values() -> None:
    frame = graph_frame()
    other = frame.assign(target_chembl_id="T2", assay_context_key="CTX2", molecule_id=lambda value: "X" + value.molecule_id)
    candidates = pd.concat([frame, other.iloc[:4]], ignore_index=True)
    selected_a, _, decision_a = select_contexts_and_targets(candidates, protocol())
    candidates["pchembl_value"] = list(reversed(candidates.pchembl_value.tolist()))
    selected_b, _, decision_b = select_contexts_and_targets(candidates, protocol())
    assert decision_a["retained_ids"] == decision_b["retained_ids"] == ["T1"]
    assert set(selected_a.activity_id) == set(selected_b.activity_id)


def test_partition_integrity_rejects_compound_or_publication_leakage() -> None:
    frame = graph_frame()
    frame["partition"] = "development"
    future = frame.iloc[[0]].copy()
    future["partition"] = "future_evaluation"
    future["activity_id"] = "future"
    with pytest.raises(AssertionError, match="Partition leakage"):
        assert_partition_integrity(pd.concat([frame, future], ignore_index=True))


def test_partition_integrity_accepts_disjoint_temporal_components() -> None:
    development = graph_frame()
    development["partition"] = "development"
    development["component_id"] = "development-component"
    future = graph_frame().replace({"T1": "T2", "M1": "X1", "M2": "X2", "M3": "X3", "M4": "X4", "D1": "F1", "D2": "F2", "D3": "F3"})
    future["partition"] = "future_evaluation"
    future["component_id"] = "future-component"
    future["document_year"] = 2021
    future["activity_id"] = "F" + future.activity_id
    future["assay_id"] = "F" + future.assay_id
    assert assert_partition_integrity(pd.concat([development, future], ignore_index=True)) == {
        "compound_overlap": 0, "publication_overlap": 0, "assay_overlap": 0, "activity_overlap": 0,
    }


def test_locked_ledger_reads_only_selected_partition(tmp_path: Path) -> None:
    (tmp_path / "assay_protocol_enrichment").mkdir()
    (tmp_path / "independent_assay_validation").mkdir()
    (tmp_path / "assay_disagreement_expansion" / "final").mkdir(parents=True)
    pd.DataFrame([
        {"partition": "development", "molecule_id": "DEV", "heldout_activity_ids": "1", "heldout_assay_ids": "A1", "heldout_document_ids": "D1"},
        {"partition": "evaluation", "molecule_id": "LOCKED", "heldout_activity_ids": "2;3", "heldout_assay_ids": "A2", "heldout_document_ids": "D2"},
    ]).to_csv(tmp_path / "assay_protocol_enrichment" / "enriched_assay_disagreement_events.csv", index=False)
    pd.DataFrame([{"molecule_id": "I", "activity_ids": "4", "assay_ids": "A4", "document_ids": "D4"}]).to_csv(tmp_path / "independent_assay_validation" / "independent_dataset.csv", index=False)
    pd.DataFrame([{"molecule_id": "E", "activity_ids": "5", "assay_ids": "A5", "document_ids": "D5"}]).to_csv(tmp_path / "assay_disagreement_expansion" / "final" / "enriched_expansion_events.csv", index=False)
    ledger, files = load_locked_ledger(tmp_path)
    assert "LOCKED" in ledger["molecules"] and "DEV" not in ledger["molecules"]
    assert ledger["activities"] == {"2", "3", "4", "5"}
    assert len(files) == 3
