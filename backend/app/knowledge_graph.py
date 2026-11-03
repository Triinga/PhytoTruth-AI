from __future__ import annotations

import hashlib
import json
import sqlite3
from collections.abc import Iterable, Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .structures import standardize_smiles


SCHEMA_VERSION = 1


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _iter_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {path}:{line_number}") from exc


def _load_specifications(path: Path) -> list[dict[str, Any]]:
    value = json.loads(path.read_text(encoding="utf-8"))
    targets = value.get("targets") if isinstance(value, dict) else value
    if not isinstance(targets, list) or not targets:
        raise ValueError("Knowledge-graph specifications require a non-empty targets list")
    output: list[dict[str, Any]] = []
    for index, item in enumerate(targets):
        if not isinstance(item, dict):
            raise ValueError(f"Target specification {index} must be an object")
        target_id = item.get("target_chembl_id") or item.get("target_id")
        activities = item.get("activities")
        if not target_id or not activities:
            raise ValueError(f"Target specification {index} requires target_chembl_id and activities")
        activity_path = Path(activities)
        if not activity_path.is_absolute():
            activity_path = (path.parent / activity_path).resolve()
        diseases = item.get("diseases", [])
        if isinstance(diseases, (str, dict)):
            diseases = [diseases]
        output.append({**item, "target_chembl_id": str(target_id), "activities": activity_path, "diseases": diseases})
    return output


class KnowledgeGraph:
    """Portable semantic graph stored in SQLite with evidence on every edge."""

    def __init__(self, path: Path, *, read_only: bool = False):
        self.path = Path(path)
        if read_only:
            self.connection = sqlite3.connect(f"file:{self.path.as_posix()}?mode=ro", uri=True)
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.connection = sqlite3.connect(self.path)
        self.connection.row_factory = sqlite3.Row
        if not read_only:
            self.connection.execute("PRAGMA journal_mode=WAL")
            self.connection.execute("PRAGMA foreign_keys=ON")
            self._create_schema()

    def _create_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS graph_metadata (
                key TEXT PRIMARY KEY,
                value_json TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS nodes (
                node_id INTEGER PRIMARY KEY,
                node_type TEXT NOT NULL,
                external_id TEXT NOT NULL,
                name TEXT,
                properties_json TEXT NOT NULL DEFAULT '{}',
                UNIQUE(node_type, external_id)
            );
            CREATE TABLE IF NOT EXISTS edges (
                edge_id INTEGER PRIMARY KEY,
                source_node_id INTEGER NOT NULL REFERENCES nodes(node_id),
                predicate TEXT NOT NULL,
                target_node_id INTEGER NOT NULL REFERENCES nodes(node_id),
                evidence_source TEXT NOT NULL,
                evidence_id TEXT NOT NULL,
                properties_json TEXT NOT NULL DEFAULT '{}',
                UNIQUE(source_node_id, predicate, target_node_id, evidence_source, evidence_id)
            );
            CREATE INDEX IF NOT EXISTS idx_nodes_type_id ON nodes(node_type, external_id);
            CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source_node_id, predicate);
            CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target_node_id, predicate);
            """
        )
        self.connection.execute(
            "INSERT OR REPLACE INTO graph_metadata VALUES (?, ?)",
            ("schema_version", json.dumps(SCHEMA_VERSION)),
        )
        self.connection.commit()

    def close(self) -> None:
        self.connection.commit()
        self.connection.close()

    def __enter__(self) -> "KnowledgeGraph":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def add_node(
        self,
        node_type: str,
        external_id: str,
        *,
        name: str | None = None,
        properties: dict[str, Any] | None = None,
    ) -> int:
        properties_json = json.dumps(properties or {}, sort_keys=True, ensure_ascii=False)
        self.connection.execute(
            """INSERT INTO nodes(node_type, external_id, name, properties_json)
               VALUES (?, ?, ?, ?)
               ON CONFLICT(node_type, external_id) DO UPDATE SET
                 name=COALESCE(excluded.name, nodes.name),
                 properties_json=CASE WHEN excluded.properties_json='{}'
                   THEN nodes.properties_json ELSE excluded.properties_json END""",
            (node_type, external_id, name, properties_json),
        )
        return int(self.connection.execute(
            "SELECT node_id FROM nodes WHERE node_type=? AND external_id=?",
            (node_type, external_id),
        ).fetchone()[0])

    def add_edge(
        self,
        source_node_id: int,
        predicate: str,
        target_node_id: int,
        *,
        evidence_source: str,
        evidence_id: str,
        properties: dict[str, Any] | None = None,
    ) -> None:
        self.connection.execute(
            """INSERT OR IGNORE INTO edges
               (source_node_id, predicate, target_node_id, evidence_source, evidence_id, properties_json)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                source_node_id,
                predicate,
                target_node_id,
                evidence_source,
                evidence_id,
                json.dumps(properties or {}, sort_keys=True, ensure_ascii=False),
            ),
        )

    def neighborhood(self, node_type: str, external_id: str, *, limit: int = 200) -> dict[str, Any] | None:
        row = self.connection.execute(
            "SELECT * FROM nodes WHERE node_type=? AND external_id=?", (node_type, external_id)
        ).fetchone()
        if row is None:
            return None
        edges = self.connection.execute(
            """SELECT e.edge_id, e.predicate, e.evidence_source, e.evidence_id, e.properties_json,
                      s.node_type source_type, s.external_id source_external_id, s.name source_name,
                      t.node_type target_type, t.external_id target_external_id, t.name target_name
               FROM edges e JOIN nodes s ON s.node_id=e.source_node_id
               JOIN nodes t ON t.node_id=e.target_node_id
               WHERE e.source_node_id=? OR e.target_node_id=? ORDER BY e.edge_id LIMIT ?""",
            (row["node_id"], row["node_id"], max(1, min(limit, 1000))),
        ).fetchall()
        return {
            "node": {**dict(row), "properties": json.loads(row["properties_json"])},
            "edges": [
                {
                    **{key: value for key, value in dict(edge).items() if key != "properties_json"},
                    "properties": json.loads(edge["properties_json"]),
                }
                for edge in edges
            ],
        }

    def summary(self) -> dict[str, Any]:
        node_counts = {row[0]: row[1] for row in self.connection.execute(
            "SELECT node_type, COUNT(*) FROM nodes GROUP BY node_type ORDER BY node_type"
        )}
        edge_counts = {row[0]: row[1] for row in self.connection.execute(
            "SELECT predicate, COUNT(*) FROM edges GROUP BY predicate ORDER BY predicate"
        )}
        return {
            "schema_version": SCHEMA_VERSION,
            "nodes": sum(node_counts.values()),
            "edges": sum(edge_counts.values()),
            "node_types": node_counts,
            "predicates": edge_counts,
        }


def _disease_value(value: Any) -> tuple[str, str]:
    if isinstance(value, str):
        return value, value
    if isinstance(value, dict):
        disease_id = value.get("id") or value.get("mondo_id") or value.get("name")
        if disease_id:
            return str(disease_id), str(value.get("name") or disease_id)
    raise ValueError(f"Invalid disease specification: {value!r}")


def build_knowledge_graph(specifications: Path, natural_product_database: Path, output: Path) -> dict[str, Any]:
    """Build plant-compound-target-disease-assay-publication relationships with provenance."""
    specifications = Path(specifications)
    natural_product_database = Path(natural_product_database)
    targets = _load_specifications(specifications)
    if not natural_product_database.is_file():
        raise FileNotFoundError(natural_product_database)
    if output.exists():
        output.unlink()

    np_connection = sqlite3.connect(natural_product_database)
    np_connection.row_factory = sqlite3.Row
    with KnowledgeGraph(output) as graph:
        for target in targets:
            target_id = target["target_chembl_id"]
            target_node = graph.add_node(
                "target",
                target_id,
                name=target.get("target_name") or target_id,
                properties={"organism": target.get("target_organism"), "target_type": target.get("target_type")},
            )
            for disease in target["diseases"]:
                disease_id, disease_name = _disease_value(disease)
                disease_node = graph.add_node("disease", disease_id, name=disease_name)
                graph.add_edge(
                    target_node,
                    "ASSOCIATED_WITH_DISEASE",
                    disease_node,
                    evidence_source="target-specification",
                    evidence_id=f"{target_id}:{disease_id}",
                    properties={"relationship": target.get("disease_relationship", "research_target")},
                )

            activity_path = target["activities"]
            if not activity_path.is_file():
                raise FileNotFoundError(activity_path)
            for record in _iter_jsonl(activity_path):
                smiles = record.get("canonical_smiles")
                if not smiles:
                    continue
                try:
                    structure = standardize_smiles(smiles)
                except (ValueError, RuntimeError):
                    continue
                compound_node = graph.add_node(
                    "compound",
                    structure.standard_inchi_key,
                    name=record.get("molecule_chembl_id") or structure.standard_inchi_key,
                    properties={
                        "canonical_smiles": structure.canonical_smiles,
                        "standard_inchi": structure.standard_inchi,
                        "molecule_chembl_id": record.get("molecule_chembl_id"),
                    },
                )
                activity_id = str(record.get("source_record_id") or record.get("activity_id") or "unidentified")
                graph.add_edge(
                    compound_node,
                    "MEASURED_AGAINST",
                    target_node,
                    evidence_source="ChEMBL",
                    evidence_id=activity_id,
                    properties={
                        "pactivity": record.get("computed_pactivity"),
                        "standard_type": record.get("standard_type"),
                        "standard_relation": record.get("standard_relation"),
                        "assay_context_key": record.get("assay_context_key"),
                    },
                )
                assay_id = record.get("assay_chembl_id")
                if assay_id:
                    assay_node = graph.add_node("assay", str(assay_id), name=str(assay_id))
                    graph.add_edge(compound_node, "TESTED_IN", assay_node, evidence_source="ChEMBL", evidence_id=activity_id)
                    graph.add_edge(assay_node, "MEASURES", target_node, evidence_source="ChEMBL", evidence_id=str(assay_id))
                document_id = record.get("document_chembl_id")
                if document_id:
                    document_node = graph.add_node("publication", str(document_id), name=str(document_id))
                    if assay_id:
                        graph.add_edge(document_node, "REPORTS", assay_node, evidence_source="ChEMBL", evidence_id=activity_id)

                provenance_rows = np_connection.execute(
                    """SELECT p.source, p.dataset_version, p.source_id, p.name, p.organism,
                              p.doi, p.collection_name
                       FROM provenance p JOIN compounds c USING(compound_id)
                       WHERE c.standard_inchi_key=? ORDER BY p.source, p.source_id""",
                    (structure.standard_inchi_key,),
                ).fetchall()
                for provenance in provenance_rows:
                    source_node = graph.add_node(
                        "natural_product_record",
                        f"{provenance['source']}:{provenance['source_id']}",
                        name=provenance["name"] or provenance["source_id"],
                        properties={
                            "database": provenance["source"],
                            "dataset_version": provenance["dataset_version"],
                            "doi": provenance["doi"],
                            "collection": provenance["collection_name"],
                        },
                    )
                    graph.add_edge(
                        compound_node,
                        "VERIFIED_BY_NATURAL_PRODUCT_SOURCE",
                        source_node,
                        evidence_source=provenance["source"],
                        evidence_id=provenance["source_id"],
                    )
                    if provenance["organism"]:
                        plant_node = graph.add_node(
                            "organism",
                            provenance["organism"],
                            name=provenance["organism"],
                            properties={"role": "reported_natural_source", "plant_status": "not_taxonomically_asserted"},
                        )
                        graph.add_edge(
                            plant_node,
                            "REPORTED_SOURCE_OF",
                            compound_node,
                            evidence_source=provenance["source"],
                            evidence_id=provenance["source_id"],
                            properties={"dataset_version": provenance["dataset_version"]},
                        )

        metadata = {
            "created_at": datetime.now(UTC).isoformat(),
            "specifications": str(specifications),
            "specifications_sha256": _sha256(specifications),
            "natural_product_database": str(natural_product_database),
            "natural_product_database_sha256": _sha256(natural_product_database),
            "target_count": len(targets),
            "semantics": "research evidence graph; edges do not imply clinical efficacy, safety, or causality",
        }
        graph.connection.execute(
            "INSERT OR REPLACE INTO graph_metadata VALUES (?, ?)", ("build", json.dumps(metadata, sort_keys=True))
        )
        graph.connection.commit()
        summary = {**graph.summary(), **metadata, "output": str(output)}
    np_connection.close()
    output.with_suffix(output.suffix + ".manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary

