from __future__ import annotations

import json
from pathlib import Path

from app.ingest import IngestionStats, NaturalProductIndex
from app.knowledge_graph import KnowledgeGraph, build_knowledge_graph
from app.screening import screen_natural_products


def _natural_product_database(path: Path) -> Path:
    stats = IngestionStats()
    records = [
        {
            "source": "NPASS", "dataset_version": "fixture", "source_file": "fixture.csv",
            "source_line": 2, "source_id": "NPC-1", "smiles": "CCOc1ccc(O)cc1",
            "name": "fixture phenol", "organism": "Salvia testensis", "doi": "10.1/test",
            "collection": "fixture", "parse_error": None, "raw": {"fixture": 1},
        },
        {
            "source": "COCONUT", "dataset_version": "fixture", "source_file": "fixture.csv",
            "source_line": 3, "source_id": "CNP-2", "smiles": "CCN(CC)CC",
            "name": "fixture amine", "organism": None, "doi": None,
            "collection": "fixture", "parse_error": None, "raw": {"fixture": 2},
        },
    ]
    with NaturalProductIndex(path) as index:
        for record in records:
            index.add(record, stats)
    return path


def test_builds_semantic_evidence_graph(tmp_path: Path) -> None:
    database = _natural_product_database(tmp_path / "natural.sqlite")
    activities = tmp_path / "activities.jsonl"
    activities.write_text(json.dumps({
        "canonical_smiles": "CCOc1ccc(O)cc1", "molecule_chembl_id": "CHEMBL_FIXTURE",
        "source_record_id": "ACT-1", "computed_pactivity": 7.1, "standard_type": "IC50",
        "standard_relation": "=", "assay_context_key": "ctx", "assay_chembl_id": "ASSAY-1",
        "document_chembl_id": "DOC-1",
    }) + "\n", encoding="utf-8")
    specs = tmp_path / "graph-specs.json"
    specs.write_text(json.dumps({"targets": [{
        "target_chembl_id": "CHEMBL335", "target_name": "PTP1B", "activities": str(activities),
        "diseases": [{"id": "MONDO:0005148", "name": "type 2 diabetes mellitus"}],
    }]}), encoding="utf-8")

    output = tmp_path / "knowledge.sqlite"
    summary = build_knowledge_graph(specs, database, output)
    assert summary["node_types"]["compound"] == 1
    assert summary["node_types"]["organism"] == 1
    assert summary["predicates"]["MEASURED_AGAINST"] == 1
    assert summary["predicates"]["REPORTED_SOURCE_OF"] == 1

    with KnowledgeGraph(output, read_only=True) as graph:
        neighborhood = graph.neighborhood("target", "CHEMBL335")
    assert neighborhood is not None
    assert {edge["predicate"] for edge in neighborhood["edges"]} >= {
        "ASSOCIATED_WITH_DISEASE", "MEASURED_AGAINST", "MEASURES"
    }


class _FakePredictor:
    def __init__(self, value: float):
        self.value = value

    def predict(self, smiles: str) -> dict:
        return {
            "canonical_smiles": smiles,
            "predicted_pchembl": self.value + (0.2 if "c1" in smiles else 0.0),
            "prediction_interval": {"lower": self.value - 0.5, "upper": self.value + 0.5, "confidence_level": 0.9, "method": "fixture"},
            "applicability_domain": {"inside": True, "status": "IN_DOMAIN_RESEARCH_ONLY", "nearest_training_tanimoto": 0.8, "threshold": 0.5, "warning": None},
            "descriptors": {},
            "model": {"name": "fixture", "version": "1", "format": "fixture", "target": "fixture", "training_structures": 10, "training_scaffolds": 5},
            "kind": "predicted",
            "limitations": ["fixture"],
        }


class _MixedDomainPredictor(_FakePredictor):
    def predict(self, smiles: str) -> dict:
        result = super().predict(smiles)
        inside = "c1" in smiles
        result["applicability_domain"] = {
            "inside": inside,
            "status": "IN_DOMAIN_RESEARCH_ONLY" if inside else "OUT_OF_DOMAIN_RESEARCH_ONLY",
            "nearest_training_tanimoto": 0.8 if inside else 0.2,
            "threshold": 0.5,
            "warning": None if inside else "Low similarity to training chemistry",
        }
        return result


def test_streaming_multitarget_screening_and_graph_enrichment(tmp_path: Path) -> None:
    database = _natural_product_database(tmp_path / "natural.sqlite")
    specs = tmp_path / "screen-specs.json"
    specs.write_text(json.dumps({
        "disease": {"id": "MONDO:TEST", "name": "fixture disease"},
        "targets": [
            {"target_chembl_id": "CHEMBL_A", "model_dir": ".", "weight": 2},
            {"target_chembl_id": "CHEMBL_B", "model_dir": ".", "weight": 1},
        ],
        "filters": {"min_molecular_weight": 0},
    }), encoding="utf-8")
    graph_path = tmp_path / "knowledge.sqlite"
    with KnowledgeGraph(graph_path):
        pass
    output = tmp_path / "screening"

    summary = screen_natural_products(
        specs, database, output, knowledge_graph=graph_path, top_k=2,
        predictors={"CHEMBL_A": _FakePredictor(7.0), "CHEMBL_B": _FakePredictor(6.5)},
    )
    assert summary["screened_compounds"] == 2
    assert summary["ranked_candidates"] == 2
    rows = [json.loads(line) for line in (output / "ranked_candidates.jsonl").read_text(encoding="utf-8").splitlines()]
    assert [row["rank"] for row in rows] == [1, 2]
    assert len(rows[0]["target_predictions"]) == 2
    assert rows[0]["research_status"] == "PRIORITIZED_FOR_EXPERIMENTAL_VALIDATION"

    with KnowledgeGraph(graph_path, read_only=True) as graph:
        disease = graph.neighborhood("disease", "MONDO:TEST")
    assert disease is not None
    assert sum(edge["predicate"] == "PRIORITIZED_FOR_DISEASE" for edge in disease["edges"]) == 2


def test_screening_can_reserve_an_out_of_domain_review_cohort(tmp_path: Path) -> None:
    database = _natural_product_database(tmp_path / "natural.sqlite")
    specs = tmp_path / "screen-specs.json"
    specs.write_text(json.dumps({
        "disease": {"id": "MONDO:TEST", "name": "fixture disease"},
        "targets": [{"target_chembl_id": "CHEMBL_A", "model_dir": ".", "weight": 1}],
        "filters": {"min_molecular_weight": 0},
    }), encoding="utf-8")

    summary = screen_natural_products(
        specs, database, tmp_path / "screening", top_k=2, out_of_domain_quota=1,
        predictors={"CHEMBL_A": _MixedDomainPredictor(7.0)},
    )

    rows = [json.loads(line) for line in (tmp_path / "screening" / "ranked_candidates.jsonl").read_text(encoding="utf-8").splitlines()]
    assert summary["in_domain_candidates"] == 1
    assert summary["out_of_domain_candidates"] == 1
    assert {row["research_status"] for row in rows} == {
        "PRIORITIZED_FOR_EXPERIMENTAL_VALIDATION", "OUT_OF_DOMAIN_REVIEW_REQUIRED"
    }

