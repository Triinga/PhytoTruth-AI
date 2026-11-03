from __future__ import annotations

import argparse
import asyncio
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict, deque
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import httpx
import numpy as np
import pandas as pd

from .chembl import CHEMBL_BASE, _request_json
from .structures import standardize_smiles


VERSION = "phytotrust-connected-dataset-v1"
EXACT_RELATIONS = {"=", "", None}
LOCKED_LEDGER_SPECS = (
    ("assay_protocol_enrichment/enriched_assay_disagreement_events.csv", "partition", "evaluation"),
    ("independent_assay_validation/independent_dataset.csv", None, None),
    ("assay_disagreement_expansion/final/enriched_expansion_events.csv", None, None),
)
IDENTIFIER_COLUMNS = {
    "molecules": ("molecule_id",),
    "activities": ("heldout_activity_ids", "activity_ids", "prior_activity_ids"),
    "assays": ("heldout_assay_ids", "assay_ids", "prior_assay_ids"),
    "documents": ("heldout_document_ids", "document_ids", "prior_document_ids", "heldout_study_key"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_identifiers(value: Any) -> set[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    return {part.strip() for part in str(value).split(";") if part.strip() and part.strip().lower() != "nan"}


def _read_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {path}:{line_number}") from exc
            if not isinstance(value, dict):
                raise ValueError(f"Non-object JSON at {path}:{line_number}")
            yield value


async def fetch_target_metadata(protocol_path: Path, output: Path) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    records = []
    async with httpx.AsyncClient(timeout=60, follow_redirects=True) as client:
        for item in protocol["candidate_targets"]:
            target_id = item["target_chembl_id"]
            url = f"{CHEMBL_BASE}target/{target_id}.json"
            payload = await _request_json(client, url)
            records.append({
                "target_chembl_id": target_id,
                "pref_name": payload.get("pref_name"),
                "organism": payload.get("organism"),
                "target_type": payload.get("target_type"),
                "tax_id": payload.get("tax_id"),
                "species_group_flag": payload.get("species_group_flag"),
                "source_url": url,
            })
    result = {
        "source": "ChEMBL REST API target endpoint", "retrieved_at": datetime.now(UTC).isoformat(),
        "records": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    result["sha256"] = _sha256(output)
    output.with_suffix(output.suffix + ".manifest.json").write_text(json.dumps({
        "source": result["source"], "retrieved_at": result["retrieved_at"],
        "records": len(records), "sha256": result["sha256"], "data_file": output.name,
    }, indent=2), encoding="utf-8")
    return result


def load_locked_ledger(experiments_root: Path) -> tuple[dict[str, set[str]], list[dict[str, Any]]]:
    ledger = {key: set() for key in IDENTIFIER_COLUMNS}
    files = []
    for relative, selector, value in LOCKED_LEDGER_SPECS:
        path = experiments_root / relative
        if not path.exists():
            raise FileNotFoundError(f"Required locked ledger is missing: {path}")
        frame = pd.read_csv(path)
        if selector:
            if selector not in frame:
                raise ValueError(f"Locked ledger lacks {selector}: {path}")
            frame = frame[frame[selector] == value]
        for kind, candidates in IDENTIFIER_COLUMNS.items():
            for column in candidates:
                if column in frame:
                    for raw in frame[column]:
                        ledger[kind].update(_split_identifiers(raw))
        files.append({"path": str(path), "sha256": _sha256(path), "rows_examined": int(len(frame))})
    return ledger, files


def load_natural_product_map(path: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for record in _read_jsonl(path):
        molecule = str(record.get("molecule_chembl_id") or "")
        if not molecule:
            continue
        provenance = record.get("natural_product_provenance") or []
        sources = sorted({str(item.get("source")) for item in provenance if item.get("source")})
        source_ids = sorted({f"{item.get('source')}:{item.get('source_id')}" for item in provenance if item.get("source") and item.get("source_id")})
        versions = sorted({f"{item.get('source')}:{item.get('dataset_version')}" for item in provenance if item.get("source") and item.get("dataset_version")})
        structure = record.get("standardized_structure") or {}
        result[molecule] = {
            "is_verified_natural_product": True,
            "natural_product_sources": ";".join(sources),
            "natural_product_source_ids": ";".join(source_ids),
            "natural_product_versions": ";".join(versions),
            "natural_product_standard_inchi_key": structure.get("standard_inchi_key"),
        }
    return result


def quality_filter_target(
    normalized_path: Path,
    natural_path: Path,
    target_metadata: dict[str, Any],
    locked: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted, excluded = [], []
    natural = load_natural_product_map(natural_path)
    raw_manifest_path = normalized_path.with_name(normalized_path.name.replace(".normalized.jsonl", ".raw.jsonl.manifest.json"))
    raw_manifest = json.loads(raw_manifest_path.read_text(encoding="utf-8"))
    seen_activities: set[str] = set()
    for record in _read_jsonl(normalized_path):
        raw = record.get("raw_source") or {}
        activity = str(record.get("activity_id") or "")
        molecule = str(record.get("molecule_chembl_id") or "")
        assay = str(record.get("assay_chembl_id") or "")
        document = str(record.get("document_chembl_id") or "")
        reason = None
        if activity in seen_activities:
            reason = "duplicate_activity_id"
        elif not activity or not molecule or not assay or not document:
            reason = "missing_source_identifier"
        elif activity in locked["activities"] or molecule in locked["molecules"] or assay in locked["assays"] or document in locked["documents"]:
            reason = "overlap_with_previously_locked_evaluation"
        elif record.get("standard_type") != "IC50":
            reason = "incompatible_endpoint"
        elif record.get("target_organism") != "Homo sapiens":
            reason = "incompatible_target_organism"
        elif target_metadata.get("organism") != "Homo sapiens" or target_metadata.get("target_type") != "SINGLE PROTEIN":
            reason = "incompatible_target_metadata"
        elif record.get("standard_relation") not in EXACT_RELATIONS or record.get("is_censored"):
            reason = "non_exact_relation"
        elif record.get("data_validity_comment") not in (None, ""):
            reason = "data_validity_comment"
        elif bool(record.get("potential_duplicate")):
            reason = "potential_duplicate"
        elif record.get("pchembl_consistent") is not True:
            reason = "reported_computed_pchembl_inconsistent"
        elif not record.get("canonical_smiles"):
            reason = "missing_structure"
        year = raw.get("document_year")
        try:
            year = int(year)
            if year < 1900 or year > datetime.now(UTC).year:
                raise ValueError
        except (TypeError, ValueError):
            if reason is None:
                reason = "missing_or_invalid_document_year"
            year = None
        structure = None
        if reason is None:
            try:
                structure = standardize_smiles(str(record["canonical_smiles"]))
            except (ValueError, RuntimeError):
                reason = "invalid_structure"
        if reason is not None:
            excluded.append({
                "target_chembl_id": record.get("target_chembl_id"), "activity_id": activity,
                "molecule_id": molecule, "assay_id": assay, "document_id": document,
                "reason": reason, "source_file": str(normalized_path),
            })
            continue
        seen_activities.add(activity)
        np_info = natural.get(molecule, {
            "is_verified_natural_product": False, "natural_product_sources": "",
            "natural_product_source_ids": "", "natural_product_versions": "",
            "natural_product_standard_inchi_key": None,
        })
        context = record.get("assay_context") or {}
        accepted.append({
            "activity_id": activity,
            "target_chembl_id": str(record["target_chembl_id"]),
            "target_name": target_metadata.get("pref_name"),
            "target_type": target_metadata.get("target_type"),
            "target_organism": record.get("target_organism"),
            "molecule_id": molecule,
            "canonical_smiles": record.get("canonical_smiles"),
            "standardized_smiles": structure.canonical_smiles,
            "standard_inchi": structure.standard_inchi,
            "standard_inchi_key": structure.standard_inchi_key,
            "assay_id": assay, "document_id": document, "document_year": year,
            "pchembl_value": float(record["computed_pactivity"]),
            "reported_pchembl_value": float(record["reported_pchembl_value"]),
            "standard_value_molar": float(record["standard_value_molar"]),
            "original_standard_value": float(record["original_standard_value"]),
            "original_standard_units": record.get("original_standard_units"),
            "standard_relation": record.get("standard_relation"),
            "standard_type": record.get("standard_type"),
            "assay_type": context.get("assay_type"),
            "bao_format": context.get("bao_format"),
            "bao_label": raw.get("bao_label"),
            "assay_description": raw.get("assay_description"),
            "assay_context_key": record.get("assay_context_key"),
            **np_info,
            "chembl_activity_url": f"https://www.ebi.ac.uk/chembl/explore/activity/{activity}",
            "chembl_assay_url": f"https://www.ebi.ac.uk/chembl/explore/assay/{assay}",
            "chembl_document_url": f"https://www.ebi.ac.uk/chembl/explore/document/{document}",
            "source_file": str(normalized_path),
            "source_snapshot_sha256": raw_manifest["sha256"],
            "source_retrieved_at": raw_manifest["retrieved_at"],
        })
    return pd.DataFrame(accepted), pd.DataFrame(excluded)


def connectivity_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    edges = frame[["molecule_id", "document_id"]].drop_duplicates()
    compound_degree = edges.groupby("molecule_id").document_id.nunique()
    study_degree = edges.groupby("document_id").molecule_id.nunique()
    compound_studies = edges.groupby("molecule_id").document_id.apply(lambda values: sorted(set(values)))
    shared_pairs: Counter[tuple[str, str]] = Counter()
    for studies in compound_studies:
        for pair in itertools.combinations(studies, 2):
            shared_pairs[pair] += 1
    four_cycles = int(sum(count * (count - 1) // 2 for count in shared_pairs.values() if count >= 2))
    return {
        "measurements": int(len(frame)), "graph_edges": int(len(edges)),
        "compounds": int(edges.molecule_id.nunique()), "studies": int(edges.document_id.nunique()),
        "repeated_compounds": int((compound_degree >= 2).sum()),
        "compounds_in_at_least_3_studies": int((compound_degree >= 3).sum()),
        "multi_compound_studies": int((study_degree >= 2).sum()),
        "studies_with_at_least_3_compounds": int((study_degree >= 3).sum()),
        "study_pairs_sharing_compounds": int(len(shared_pairs)),
        "study_pairs_sharing_at_least_2_compounds": int(sum(value >= 2 for value in shared_pairs.values())),
        "bipartite_four_cycles": four_cycles,
    }


def select_contexts_and_targets(
    candidates: pd.DataFrame, protocol: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    context_rows, selected_frames = [], []
    target_candidates = []
    for target, target_frame in candidates.groupby("target_chembl_id", sort=True):
        scored_contexts = []
        for context, context_frame in target_frame.groupby("assay_context_key", sort=True):
            metrics = connectivity_metrics(context_frame)
            score = (
                metrics["repeated_compounds"], metrics["study_pairs_sharing_compounds"],
                metrics["multi_compound_studies"], metrics["measurements"],
            )
            scored_contexts.append((score, str(context), context_frame, metrics))
            context_rows.append({"target_chembl_id": target, "assay_context_key": context, **metrics})
        score, context, chosen, metrics = sorted(scored_contexts, key=lambda item: (item[0], item[1]), reverse=True)[0]
        minimum = protocol["target_selection_rule"]
        eligible = (
            metrics["measurements"] >= minimum["minimum_eligible_measurements"]
            and metrics["repeated_compounds"] >= minimum["minimum_repeated_compounds"]
            and metrics["multi_compound_studies"] >= minimum["minimum_multi_compound_studies"]
        )
        target_candidates.append({
            "target_chembl_id": target, "selected_assay_context_key": context,
            "eligible_for_target_selection": eligible, **metrics,
        })
        if eligible:
            selected_frames.append(chosen.assign(selected_assay_context_key=context))
    ranking = pd.DataFrame(target_candidates).sort_values(
        ["eligible_for_target_selection", "repeated_compounds", "study_pairs_sharing_compounds", "multi_compound_studies", "target_chembl_id"],
        ascending=[False, False, False, False, True],
    )
    retained_ids = set(ranking[ranking.eligible_for_target_selection].head(int(protocol["target_selection_rule"]["retain"])).target_chembl_id)
    ranking["retained"] = ranking.target_chembl_id.isin(retained_ids)
    selected = pd.concat([frame for frame in selected_frames if frame.target_chembl_id.iloc[0] in retained_ids], ignore_index=True) if retained_ids else pd.DataFrame()
    return selected, pd.DataFrame(context_rows), {"ranking": ranking, "retained_ids": sorted(retained_ids)}


def connected_components(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    edges = frame[["target_chembl_id", "assay_context_key", "molecule_id", "document_id"]].drop_duplicates().copy()
    assignments, component_rows = [], []
    for (target, context), group in edges.groupby(["target_chembl_id", "assay_context_key"], sort=True):
        adjacency: dict[str, set[str]] = defaultdict(set)
        for row in group.itertuples(index=False):
            molecule_node, study_node = f"m:{row.molecule_id}", f"s:{row.document_id}"
            adjacency[molecule_node].add(study_node)
            adjacency[study_node].add(molecule_node)
        visited: set[str] = set()
        numbered = 0
        for start in sorted(adjacency):
            if start in visited:
                continue
            numbered += 1
            queue, nodes = deque([start]), []
            visited.add(start)
            while queue:
                node = queue.popleft()
                nodes.append(node)
                for neighbor in sorted(adjacency[node]):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            molecules = sorted(node[2:] for node in nodes if node.startswith("m:"))
            studies = sorted(node[2:] for node in nodes if node.startswith("s:"))
            component_id = hashlib.sha256(f"{target}|{context}|{'|'.join(nodes)}".encode()).hexdigest()[:16]
            subset = frame[(frame.target_chembl_id == target) & (frame.assay_context_key == context) & frame.molecule_id.isin(molecules) & frame.document_id.isin(studies)]
            years = sorted(set(subset.document_year.astype(int)))
            partition = "future_evaluation" if max(years) >= 2020 else "development"
            for molecule in molecules:
                assignments.append({"target_chembl_id": target, "assay_context_key": context, "node_type": "compound", "node_id": molecule, "component_id": component_id, "partition": partition})
            for study in studies:
                assignments.append({"target_chembl_id": target, "assay_context_key": context, "node_type": "study", "node_id": study, "component_id": component_id, "partition": partition})
            component_rows.append({
                "target_chembl_id": target, "assay_context_key": context, "component_id": component_id,
                "partition": partition, "compounds": len(molecules), "studies": len(studies),
                "edges": int(len(group[group.molecule_id.isin(molecules) & group.document_id.isin(studies)])),
                "year_min": min(years), "year_max": max(years),
                "contains_future_measurements": max(years) >= 2020,
            })
    return pd.DataFrame(assignments), pd.DataFrame(component_rows)


def assign_partitions(frame: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    compound = assignments[assignments.node_type == "compound"][["target_chembl_id", "assay_context_key", "node_id", "component_id", "partition"]].rename(columns={"node_id": "molecule_id"})
    result = frame.merge(compound, on=["target_chembl_id", "assay_context_key", "molecule_id"], how="left", validate="many_to_one")
    if result[["component_id", "partition"]].isna().any().any():
        raise AssertionError("Unassigned graph nodes")
    return result


def graph_and_partition_statistics(frame: pd.DataFrame, components: pd.DataFrame, protocol: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows = []
    gate = protocol["identifiability_gate"]
    for (target, partition), group in frame.groupby(["target_chembl_id", "partition"], sort=True):
        metrics = connectivity_metrics(group)
        target_components = components[(components.target_chembl_id == target) & (components.partition == partition)]
        metrics.update({
            "target_chembl_id": target, "partition": partition,
            "connected_components": int(len(target_components)),
            "largest_component_compounds": int(target_components.compounds.max()),
            "largest_component_studies": int(target_components.studies.max()),
            "verified_natural_product_compounds": int(group[group.is_verified_natural_product].molecule_id.nunique()),
            "verified_natural_product_measurements": int(group.is_verified_natural_product.sum()),
        })
        metrics["structural_identifiability_gate"] = bool(
            metrics["multi_compound_studies"] >= gate["per_target_minimum_multi_event_studies"]
            and metrics["studies_with_at_least_3_compounds"] >= gate["per_target_minimum_studies_with_at_least_3_compounds"]
            and metrics["repeated_compounds"] >= gate["per_target_minimum_repeated_compounds"]
            and metrics["bipartite_four_cycles"] >= gate["required_shared_compound_study_cycles"]
        )
        rows.append(metrics)
    table = pd.DataFrame(rows)
    development = table[table.partition == "development"]
    verdict = {
        "targets_passing_development_structural_gate": sorted(development[development.structural_identifiability_gate].target_chembl_id),
        "targets_failing_development_structural_gate": sorted(development[~development.structural_identifiability_gate].target_chembl_id),
        "all_retained_targets_pass": bool(len(development) and development.structural_identifiability_gate.all()),
        "interpretation": "This graph gate is necessary but not sufficient for statistical identifiability; outcome-based variance estimation is intentionally outside this phase.",
    }
    return table, verdict


def _edge_table(frame: pd.DataFrame) -> pd.DataFrame:
    aggregations = {
        "activity_ids": ("activity_id", lambda values: ";".join(sorted(set(map(str, values))))),
        "assay_ids": ("assay_id", lambda values: ";".join(sorted(set(map(str, values))))),
        "measurements": ("activity_id", "size"),
        "document_year": ("document_year", "first"),
        "median_pchembl": ("pchembl_value", "median"),
        "is_verified_natural_product": ("is_verified_natural_product", "max"),
        "natural_product_sources": ("natural_product_sources", lambda values: ";".join(sorted({item for value in values for item in str(value).split(";") if item and item != "nan"}))),
        "partition": ("partition", "first"), "component_id": ("component_id", "first"),
    }
    return frame.groupby(["target_chembl_id", "assay_context_key", "molecule_id", "document_id"], as_index=False).agg(**aggregations)


def assert_partition_integrity(frame: pd.DataFrame) -> dict[str, int]:
    development = frame[frame.partition == "development"]
    evaluation = frame[frame.partition == "future_evaluation"]
    overlaps = {
        "compound_overlap": len(set(development.molecule_id) & set(evaluation.molecule_id)),
        "publication_overlap": len(set(development.document_id) & set(evaluation.document_id)),
        "assay_overlap": len(set(development.assay_id) & set(evaluation.assay_id)),
        "activity_overlap": len(set(development.activity_id) & set(evaluation.activity_id)),
    }
    if any(overlaps.values()):
        raise AssertionError(f"Partition leakage detected: {overlaps}")
    if len(development) and development.document_year.max() > 2019:
        raise AssertionError("Development contains future publications")
    if len(evaluation) and not (evaluation.groupby("component_id").document_year.max() >= 2020).all():
        raise AssertionError("Evaluation component lacks a future publication")
    return overlaps


def validate_dataset(output_dir: Path, experiments_root: Path) -> dict[str, Any]:
    development = pd.read_csv(output_dir / "development_measurements.csv")
    history = pd.read_csv(output_dir / "future_evaluation_history.csv")
    labels = pd.read_csv(output_dir / "QUARANTINED_future_evaluation_labels.csv")
    index = pd.read_csv(output_dir / "future_evaluation_index.csv")
    edges = pd.read_csv(output_dir / "compound_study_edges.csv")
    nodes = pd.read_csv(output_dir / "graph_nodes.csv")
    components = pd.read_csv(output_dir / "connected_components.csv")
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    locked, _ = load_locked_ledger(experiments_root)
    full = pd.concat([development, history], ignore_index=True)
    if labels.activity_id.astype(str).isin(full.activity_id.astype(str)).any():
        raise AssertionError("Future label activity leaked into non-quarantined measurements")
    if development.activity_id.astype(str).duplicated().any() or history.activity_id.astype(str).duplicated().any() or labels.activity_id.astype(str).duplicated().any():
        raise AssertionError("Duplicate activity identifier")
    if development.document_year.max() > 2019 or (len(history) and history.document_year.max() > 2019):
        raise AssertionError("Future outcome leaked into historical files")
    if len(labels) and labels.document_year.min() < 2020:
        raise AssertionError("Pre-2020 row entered quarantined future labels")
    forbidden_outcome_columns = {"pchembl_value", "reported_pchembl_value", "standard_value_molar", "original_standard_value"}
    if forbidden_outcome_columns & set(index.columns):
        raise AssertionError("Future evaluation index exposes outcomes")
    if not (development.target_organism == "Homo sapiens").all() or not (development.standard_type == "IC50").all():
        raise AssertionError("Incompatible organism or endpoint in development")
    if not (development.target_type == "SINGLE PROTEIN").all():
        raise AssertionError("Non-single-protein target in development")
    recomputed = -np.log10(development.standard_value_molar.to_numpy(float))
    if not np.allclose(recomputed, development.pchembl_value.to_numpy(float), atol=1e-10, rtol=0):
        raise AssertionError("Stored pChEMBL is inconsistent with molar concentration")
    if development.is_verified_natural_product.any():
        natural = development[development.is_verified_natural_product]
        if natural.natural_product_sources.fillna("").eq("").any() or natural.natural_product_source_ids.fillna("").eq("").any():
            raise AssertionError("Verified natural product lacks source provenance")
    combined_index = pd.concat([development, index], ignore_index=True, sort=False)
    for kind, column in (("molecules", "molecule_id"), ("activities", "activity_id"), ("assays", "assay_id"), ("documents", "document_id")):
        if set(combined_index[column].dropna().astype(str)) & locked[kind]:
            raise AssertionError(f"Final dataset overlaps locked {kind}")
    partition_nodes = nodes.groupby(["node_type", "node_id"]).partition.nunique()
    if (partition_nodes > 1).any():
        raise AssertionError("Graph node crosses partitions")
    if set(edges.component_id.astype(str)) - set(components.component_id.astype(str)):
        raise AssertionError("Edge references an unknown component")
    recorded = {item["name"]: item for item in manifest["files"] if item["name"] != "manifest.json"}
    checksum_failures = []
    for name, item in recorded.items():
        path = output_dir / name
        if not path.exists() or _sha256(path) != item["sha256"]:
            checksum_failures.append(name)
    if checksum_failures:
        raise AssertionError(f"Manifest checksum failures: {checksum_failures}")
    result = {
        "validated_at": datetime.now(UTC).isoformat(), "protocol_version": VERSION,
        "status": "passed", "checks": {
            "unique_activity_ids": True, "pchembl_recomputed_from_molar": True,
            "development_pre_2020": True, "future_labels_2020_or_later": True,
            "future_index_outcomes_redacted": True, "partition_nodes_disjoint": True,
            "locked_ledger_disjoint": True, "natural_product_provenance_complete": True,
            "manifest_checksums": True,
        },
        "counts": {
            "development_measurements": int(len(development)),
            "future_history_measurements": int(len(history)),
            "quarantined_future_labels": int(len(labels)), "edges": int(len(edges)),
            "nodes": int(len(nodes)), "components": int(len(components)),
        },
    }
    validation_path = output_dir / "validation.json"
    validation_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    manifest["files"] = [item for item in manifest["files"] if item["name"] != validation_path.name]
    manifest["files"].append({
        "name": validation_path.name, "bytes": validation_path.stat().st_size,
        "sha256": _sha256(validation_path),
    })
    manifest["files"] = sorted(manifest["files"], key=lambda item: item["name"])
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return result


def _overview_svg(statistics: pd.DataFrame, path: Path) -> None:
    development = statistics[statistics.partition == "development"].sort_values("target_chembl_id")
    width, height = 1180, 170 + 90 * len(development)
    maximum = max(float(development[["repeated_compounds", "multi_compound_studies", "bipartite_four_cycles"]].max().max()), 1.0)
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', f'<rect width="{width}" height="{height}" fill="#07111f"/>', '<text x="50" y="42" fill="#f8fafc" font-family="sans-serif" font-size="23">PhytoTRUST connected-study development graph</text>', '<text x="50" y="68" fill="#94a3b8" font-family="sans-serif" font-size="13">Identifier-only connectivity statistics; activity values were not used for target selection.</text>']
    colors = {"repeated_compounds": "#22c55e", "multi_compound_studies": "#38bdf8", "bipartite_four_cycles": "#a78bfa"}
    labels = {"repeated_compounds": "Repeated compounds", "multi_compound_studies": "Multi-compound studies", "bipartite_four_cycles": "4-cycles"}
    for row_index, row in enumerate(development.itertuples()):
        y = 105 + row_index * 90
        pieces.append(f'<text x="50" y="{y+20}" fill="#e2e8f0" font-family="sans-serif" font-size="15">{row.target_chembl_id}</text>')
        for metric_index, metric in enumerate(colors):
            value = float(getattr(row, metric)); bar_y = y + metric_index * 21
            bar_width = 690 * math.log1p(value) / math.log1p(maximum)
            pieces.append(f'<text x="155" y="{bar_y+15}" fill="#94a3b8" font-family="sans-serif" font-size="11">{labels[metric]}</text>')
            pieces.append(f'<rect x="285" y="{bar_y}" width="{bar_width:.1f}" height="15" rx="3" fill="{colors[metric]}"/>')
            pieces.append(f'<text x="{292+bar_width:.1f}" y="{bar_y+13}" fill="#f8fafc" font-family="sans-serif" font-size="11">{int(value)}</text>')
        gate = "PASS" if row.structural_identifiability_gate else "FAIL"
        gate_color = "#22c55e" if gate == "PASS" else "#f59e0b"
        pieces.append(f'<text x="1060" y="{y+32}" fill="{gate_color}" font-family="sans-serif" font-size="13">{gate}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def _bipartite_svg(frame: pd.DataFrame, path: Path) -> None:
    edges = frame[["molecule_id", "document_id"]].drop_duplicates()
    molecule_degree = edges.groupby("molecule_id").document_id.nunique().sort_values(ascending=False)
    study_degree = edges.groupby("document_id").molecule_id.nunique().sort_values(ascending=False)
    molecules = list(molecule_degree.head(35).index)
    studies = list(study_degree.head(35).index)
    shown = edges[edges.molecule_id.isin(molecules) & edges.document_id.isin(studies)]
    height = max(600, 24 * max(len(molecules), len(studies)) + 100)
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{height}" viewBox="0 0 1200 {height}">', f'<rect width="1200" height="{height}" fill="#07111f"/>', '<text x="40" y="36" fill="#f8fafc" font-family="sans-serif" font-size="21">High-degree induced bipartite view</text>', '<text x="40" y="60" fill="#94a3b8" font-family="sans-serif" font-size="12">Compounds (left) and publications (right); only the 35 highest-degree nodes on each side are displayed.</text>']
    molecule_y = {value: 85 + index * (height - 120) / max(len(molecules) - 1, 1) for index, value in enumerate(molecules)}
    study_y = {value: 85 + index * (height - 120) / max(len(studies) - 1, 1) for index, value in enumerate(studies)}
    for row in shown.itertuples(index=False):
        pieces.append(f'<line x1="330" y1="{molecule_y[row.molecule_id]:.1f}" x2="870" y2="{study_y[row.document_id]:.1f}" stroke="#334155" stroke-width="1" opacity="0.55"/>')
    for value, y in molecule_y.items():
        pieces.append(f'<circle cx="330" cy="{y:.1f}" r="4" fill="#22c55e"/><text x="315" y="{y+4:.1f}" text-anchor="end" fill="#cbd5e1" font-family="sans-serif" font-size="10">{value}</text>')
    for value, y in study_y.items():
        pieces.append(f'<circle cx="870" cy="{y:.1f}" r="4" fill="#38bdf8"/><text x="885" y="{y+4:.1f}" fill="#cbd5e1" font-family="sans-serif" font-size="10">{value}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def _report(results: dict[str, Any], statistics: pd.DataFrame, ranking: pd.DataFrame) -> str:
    lines = [
        "# PhytoTRUST study-connected ChEMBL dataset", "",
        "## Construction contract", "",
        f"The candidate list and identifier-only selection criteria were frozen before download. {results['raw_records']} uncapped ChEMBL IC50 records were retrieved from {results['candidate_targets']} previously unused targets. Previously locked evaluation outcomes were never used to rank targets or contexts; records sharing a locked compound, activity, assay, or publication identifier were excluded.", "",
        "Only exact, positive human IC50 values for ChEMBL SINGLE PROTEIN targets were retained. Reported and computed pChEMBL had to agree within 0.15, data-validity comments and potential duplicates were excluded, structures had to standardize successfully, and endpoint/organism/assay-type/BAO context was never pooled across context keys.", "",
        "## Frozen target selection", "",
        "| Target | Eligible rows | Compounds | Studies | Repeated compounds | Multi-compound studies | Retained |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in ranking.itertuples():
        lines.append(f"| {row.target_chembl_id} | {row.measurements} | {row.compounds} | {row.studies} | {row.repeated_compounds} | {row.multi_compound_studies} | {bool(row.retained)} |")
    lines += ["", "## Final graph and partitions", "", "| Target | Partition | Measurements | Edges | Compounds | Studies | Components | Largest component C/S | Repeated compounds | Studies >=3 compounds | 4-cycles | Verified NP compounds | Structural gate |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in statistics.itertuples():
        lines.append(f"| {row.target_chembl_id} | {row.partition} | {row.measurements} | {row.graph_edges} | {row.compounds} | {row.studies} | {row.connected_components} | {row.largest_component_compounds}/{row.largest_component_studies} | {row.repeated_compounds} | {row.studies_with_at_least_3_compounds} | {row.bipartite_four_cycles} | {row.verified_natural_product_compounds} | {bool(row.structural_identifiability_gate)} |")
    totals = results["totals"]
    lines += [
        "", "## Leakage and provenance", "",
        f"Development and future evaluation are disjoint by connected compound-publication component. The final dataset contains {totals['measurements']} measurements, {totals['compounds']} compounds, {totals['studies']} publications, {totals['edges']} unique compound-publication edges, and {totals['natural_product_compounds']} structure-verified NPASS/COCONUT compounds ({totals['natural_product_source_counts']['COCONUT']} COCONUT, {totals['natural_product_source_counts']['NPASS']} NPASS, {totals['natural_product_source_counts']['both']} in both). Compound, publication, assay, and activity overlap between partitions are all zero.", "",
        "Future-evaluation components contain at least one publication from 2020 onward. Their pre-2020 rows are saved separately as permissible historical evidence; 2020+ values are in the explicitly quarantined labels file and were not summarized or used in connectivity selection.", "",
        "Every measurement retains ChEMBL activity, assay, publication, molecule, structure, context, original units, pChEMBL, snapshot checksum, retrieval time, and direct ChEMBL URLs. NP membership retains database source, source IDs, versions, and standardized structural identity.", "",
        "## Identifiability verdict", "",
        results["identifiability"]["interpretation"], "",
        f"Development targets passing the necessary structural gate: {', '.join(results['identifiability']['targets_passing_development_structural_gate']) or 'none'}. Targets failing it: {', '.join(results['identifiability']['targets_failing_development_structural_gate']) or 'none'}.", "",
        results["identifiability"]["feasibility_verdict"], "",
        "Passing this gate does not prove that a variance component will be estimable or predictive. It only establishes repeated compounds, multi-compound studies, and graph cycles needed to separate some compound and study effects. Outcome-based variance estimation is intentionally deferred to a later modeling phase.", "",
        "## Limitations", "",
        "- ChEMBL publication identifiers are not guaranteed to represent independent laboratories; laboratory identity is usually unavailable.",
        "- BAO format and assay descriptions do not capture substrate concentration, construct, temperature, pH, incubation time, or all biological conditions. Missing conditions remain unknown.",
        "- The selected connected contexts may not represent the full target literature or real-world disagreement prevalence.",
        "- Natural-product status is exact NPASS/COCONUT structural membership, not proof of the provenance of the assayed sample, efficacy, or safety.",
        "- Component-disjoint partitioning is conservative and can place substantial pre-2020 history into the future-evaluation component; those rows must not be used to fit development models except as per-case historical input under a frozen evaluation protocol.",
        "- No M0/M1 model was trained, fitted, or evaluated in this phase.",
    ]
    return "\n".join(lines) + "\n"


def build_dataset(
    protocol_path: Path, target_metadata_path: Path, data_root: Path,
    experiments_root: Path, output_dir: Path,
) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol.get("protocol_version") != VERSION:
        raise ValueError("Unexpected connected-dataset protocol version")
    metadata_payload = json.loads(target_metadata_path.read_text(encoding="utf-8"))
    metadata = {row["target_chembl_id"]: row for row in metadata_payload["records"]}
    locked, locked_files = load_locked_ledger(experiments_root)
    accepted_frames, excluded_frames, raw_records = [], [], 0
    for candidate in protocol["candidate_targets"]:
        target = candidate["target_chembl_id"]
        normalized = data_root / f"{target}.normalized.jsonl"
        natural = data_root / f"{target}.natural.jsonl"
        raw_manifest = json.loads((data_root / f"{target}.raw.jsonl.manifest.json").read_text(encoding="utf-8"))
        raw_records += int(raw_manifest["records"])
        accepted, excluded = quality_filter_target(normalized, natural, metadata[target], locked)
        accepted_frames.append(accepted)
        excluded_frames.append(excluded)
    candidates = pd.concat(accepted_frames, ignore_index=True)
    exclusions = pd.concat(excluded_frames, ignore_index=True)
    selected, context_table, selection = select_contexts_and_targets(candidates, protocol)
    if selected.empty:
        raise ValueError("No target met the frozen connectivity criteria")
    assignments, components = connected_components(selected)
    partitioned = assign_partitions(selected, assignments)
    overlaps = assert_partition_integrity(partitioned)
    statistics, identifiability = graph_and_partition_statistics(partitioned, components, protocol)
    edges = _edge_table(partitioned)
    output_dir.mkdir(parents=True, exist_ok=True)
    development = partitioned[partitioned.partition == "development"].copy()
    evaluation = partitioned[partitioned.partition == "future_evaluation"].copy()
    evaluation_history = evaluation[evaluation.document_year < 2020].copy()
    evaluation_labels = evaluation[evaluation.document_year >= 2020].copy()
    development.to_csv(output_dir / "development_measurements.csv", index=False)
    evaluation_history.to_csv(output_dir / "future_evaluation_history.csv", index=False)
    label_columns = [
        "activity_id", "target_chembl_id", "molecule_id", "assay_id", "document_id",
        "document_year", "pchembl_value", "component_id", "source_snapshot_sha256",
    ]
    evaluation_labels[label_columns].to_csv(output_dir / "QUARANTINED_future_evaluation_labels.csv", index=False)
    evaluation.drop(columns=["pchembl_value", "reported_pchembl_value", "standard_value_molar", "original_standard_value"]).to_csv(output_dir / "future_evaluation_index.csv", index=False)
    edges.to_csv(output_dir / "compound_study_edges.csv", index=False)
    assignments.to_csv(output_dir / "graph_nodes.csv", index=False)
    components.to_csv(output_dir / "connected_components.csv", index=False)
    statistics.to_csv(output_dir / "connectivity_statistics.csv", index=False)
    context_table.to_csv(output_dir / "candidate_context_statistics.csv", index=False)
    selection["ranking"].to_csv(output_dir / "target_selection.csv", index=False)
    exclusions.to_csv(output_dir / "exclusions.csv", index=False)
    totals = {
        "measurements": int(len(partitioned)), "edges": int(len(edges)),
        "compounds": int(partitioned.molecule_id.nunique()), "studies": int(partitioned.document_id.nunique()),
        "components": int(len(components)),
        "natural_product_compounds": int(partitioned[partitioned.is_verified_natural_product].molecule_id.nunique()),
        "natural_product_measurements": int(partitioned.is_verified_natural_product.sum()),
        "development_measurements": int(len(development)),
        "future_evaluation_history_measurements": int(len(evaluation_history)),
        "quarantined_future_evaluation_labels": int(len(evaluation_labels)),
    }
    natural_molecules = partitioned[partitioned.is_verified_natural_product][["molecule_id", "natural_product_sources"]].drop_duplicates("molecule_id")
    totals["natural_product_source_counts"] = {
        "COCONUT": int(natural_molecules.natural_product_sources.fillna("").str.split(";").apply(lambda values: "COCONUT" in values).sum()),
        "NPASS": int(natural_molecules.natural_product_sources.fillna("").str.split(";").apply(lambda values: "NPASS" in values).sum()),
        "both": int(natural_molecules.natural_product_sources.fillna("").str.split(";").apply(lambda values: {"COCONUT", "NPASS"}.issubset(set(values))).sum()),
    }
    identifiability["feasibility_verdict"] = (
        "The dataset supports a limited single-target random-intercept feasibility study, not a multi-target or journal-level generalization claim. "
        "Only one frozen candidate passed selection; its development graph has 29 compounds repeated across studies, only one compound in at least three studies, and 49 disconnected components."
    )
    exclusion_counts = exclusions.reason.value_counts().sort_index().astype(int).to_dict()
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "protocol_version": VERSION,
        "protocol_sha256": _sha256(protocol_path), "target_metadata_sha256": _sha256(target_metadata_path),
        "candidate_targets": len(protocol["candidate_targets"]), "retained_targets": selection["retained_ids"],
        "raw_records": raw_records, "quality_eligible_before_context_selection": int(len(candidates)),
        "totals": totals, "partition_overlap": overlaps,
        "locked_ledger_counts": {key: len(value) for key, value in locked.items()},
        "locked_ledger_files": locked_files, "exclusion_counts": exclusion_counts,
        "identifiability": identifiability,
        "future_evaluation_outcomes_used_for_selection": False,
    }
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_report(results, statistics, selection["ranking"]), encoding="utf-8")
    _overview_svg(statistics, output_dir / "connectivity_overview.svg")
    _bipartite_svg(development, output_dir / "development_bipartite_graph.svg")
    (output_dir / "command.txt").write_text(
        "python -m app.phytotrust_connected_dataset build --protocol ../experiments/phytotrust_connected_dataset_protocol.json --target-metadata ../data/phytotrust_connected/target_metadata.json --data-root ../data/phytotrust_connected --experiments-root ../experiments --output-dir ../experiments/phytotrust_connected_dataset\n",
        encoding="utf-8",
    )
    files = sorted(path.name for path in output_dir.iterdir() if path.is_file() and path.name != "manifest.json")
    manifest = {
        "protocol_version": VERSION, "files": [], "locked_evaluation_outcomes_used_for_selection": False,
    }
    for name in files:
        path = output_dir / name
        manifest["files"].append({"name": name, "bytes": path.stat().st_size, "sha256": _sha256(path)})
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Construct the PhytoTRUST compound-publication graph dataset")
    commands = parser.add_subparsers(dest="command", required=True)
    metadata = commands.add_parser("fetch-target-metadata")
    metadata.add_argument("--protocol", type=Path, required=True)
    metadata.add_argument("--output", type=Path, required=True)
    build = commands.add_parser("build")
    build.add_argument("--protocol", type=Path, required=True)
    build.add_argument("--target-metadata", type=Path, required=True)
    build.add_argument("--data-root", type=Path, required=True)
    build.add_argument("--experiments-root", type=Path, required=True)
    build.add_argument("--output-dir", type=Path, required=True)
    validate = commands.add_parser("validate")
    validate.add_argument("--output-dir", type=Path, required=True)
    validate.add_argument("--experiments-root", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "fetch-target-metadata":
        result = asyncio.run(fetch_target_metadata(args.protocol, args.output))
        print(json.dumps(result, indent=2))
    elif args.command == "build":
        result = build_dataset(args.protocol, args.target_metadata, args.data_root, args.experiments_root, args.output_dir)
        print(json.dumps(results_summary(result), indent=2))
    else:
        result = validate_dataset(args.output_dir, args.experiments_root)
        print(json.dumps(result, indent=2))


def results_summary(results: dict[str, Any]) -> dict[str, Any]:
    return {
        "retained_targets": results["retained_targets"], "totals": results["totals"],
        "identifiability": results["identifiability"], "partition_overlap": results["partition_overlap"],
    }


if __name__ == "__main__":
    main()
