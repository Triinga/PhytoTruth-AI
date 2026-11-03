from __future__ import annotations

import csv
import hashlib
import heapq
import json
import math
import sqlite3
from collections.abc import Iterator, Mapping
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Protocol

from .features import molecular_descriptors
from .knowledge_graph import KnowledgeGraph
from .prediction_service import ActivityPredictor, load_predictor


class Predictor(Protocol):
    def predict(self, smiles: str) -> dict[str, Any]: ...


DEFAULT_FILTERS: dict[str, float] = {
    "min_molecular_weight": 100.0,
    "max_molecular_weight": 650.0,
    "min_logp": -2.0,
    "max_logp": 6.0,
    "max_tpsa": 180.0,
    "max_hbd": 8.0,
    "max_hba": 14.0,
    "max_rotatable_bonds": 15.0,
    "max_abs_formal_charge": 2.0,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_specification(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or not isinstance(value.get("targets"), list) or not value["targets"]:
        raise ValueError("Screening specification requires a non-empty targets list")
    disease = value.get("disease")
    if isinstance(disease, str):
        disease = {"id": disease, "name": disease}
    if not isinstance(disease, dict) or not (disease.get("id") or disease.get("name")):
        raise ValueError("Screening specification requires disease.id or disease.name")
    targets: list[dict[str, Any]] = []
    for index, target in enumerate(value["targets"]):
        if not isinstance(target, dict):
            raise ValueError(f"Target specification {index} must be an object")
        target_id = target.get("target_chembl_id") or target.get("target_id")
        model_dir = target.get("model_dir")
        if not target_id or not model_dir:
            raise ValueError(f"Target specification {index} requires target_chembl_id and model_dir")
        resolved = Path(model_dir)
        if not resolved.is_absolute():
            resolved = (path.parent / resolved).resolve()
        weight = float(target.get("weight", 1.0))
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError(f"Target specification {index} weight must be positive")
        targets.append({**target, "target_chembl_id": str(target_id), "model_dir": resolved, "weight": weight})
    filters = {**DEFAULT_FILTERS, **value.get("filters", {})}
    return {**value, "disease": disease, "targets": targets, "filters": filters}


def _iter_compounds(
    database: Path,
    *,
    source: str | None,
    max_compounds: int | None,
) -> Iterator[dict[str, Any]]:
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    where = "WHERE p.source=?" if source else ""
    parameters: list[Any] = [source.upper()] if source else []
    limit = " LIMIT ?" if max_compounds is not None else ""
    if max_compounds is not None:
        if max_compounds < 1:
            raise ValueError("max_compounds must be positive")
        parameters.append(max_compounds)
    query = f"""
        SELECT c.standard_inchi_key, c.canonical_smiles,
               GROUP_CONCAT(DISTINCT p.source || ':' || p.source_id) provenance_ids,
               GROUP_CONCAT(DISTINCT p.organism) organisms,
               GROUP_CONCAT(DISTINCT p.name) names
        FROM compounds c JOIN provenance p USING(compound_id)
        {where}
        GROUP BY c.compound_id
        ORDER BY c.compound_id{limit}
    """
    try:
        for row in connection.execute(query, parameters):
            yield dict(row)
    finally:
        connection.close()


def _property_filter(descriptors: Mapping[str, float], filters: Mapping[str, float]) -> list[str]:
    failures: list[str] = []
    checks = (
        (descriptors["molecular_weight"] >= filters["min_molecular_weight"], "molecular_weight_below_minimum"),
        (descriptors["molecular_weight"] <= filters["max_molecular_weight"], "molecular_weight_above_maximum"),
        (descriptors["logp"] >= filters["min_logp"], "logp_below_minimum"),
        (descriptors["logp"] <= filters["max_logp"], "logp_above_maximum"),
        (descriptors["tpsa"] <= filters["max_tpsa"], "tpsa_above_maximum"),
        (descriptors["hbd"] <= filters["max_hbd"], "hbd_above_maximum"),
        (descriptors["hba"] <= filters["max_hba"], "hba_above_maximum"),
        (descriptors["rotatable_bonds"] <= filters["max_rotatable_bonds"], "rotatable_bonds_above_maximum"),
        (abs(descriptors["formal_charge"]) <= filters["max_abs_formal_charge"], "formal_charge_outside_range"),
    )
    for passed, reason in checks:
        if not passed:
            failures.append(reason)
    return failures


@lru_cache(maxsize=1)
def _alert_catalog():
    try:
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        return FilterCatalog(params)
    except (ImportError, AttributeError):
        return None


def _structural_alerts(smiles: str) -> list[str]:
    """Return PAINS/Brenk research alerts; absence is not a safety claim."""
    try:
        from rdkit import Chem

        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return ["invalid_structure"]
        catalog = _alert_catalog()
        if catalog is None:
            return []
        return sorted({entry.GetDescription() for entry in catalog.GetMatches(molecule)})
    except (ImportError, AttributeError):
        return []


def _target_score(prediction: Mapping[str, Any]) -> float:
    interval = prediction["prediction_interval"]
    width = max(0.0, float(interval["upper"]) - float(interval["lower"]))
    in_domain = bool(prediction["applicability_domain"]["inside"])
    similarity = float(prediction["applicability_domain"]["nearest_training_tanimoto"])
    return float(prediction["predicted_pchembl"]) + (0.35 if in_domain else -1.25) + 0.25 * similarity - 0.15 * width


def _candidate_row(
    compound: Mapping[str, Any],
    specification: Mapping[str, Any],
    predictors: Mapping[str, Predictor],
) -> tuple[dict[str, Any] | None, str | None]:
    smiles = str(compound["canonical_smiles"])
    descriptors = molecular_descriptors(smiles)
    filter_failures = _property_filter(descriptors, specification["filters"])
    if filter_failures:
        return None, ",".join(filter_failures)
    alerts = _structural_alerts(smiles)
    predictions: list[dict[str, Any]] = []
    weighted_score = weight_total = 0.0
    for target in specification["targets"]:
        target_id = target["target_chembl_id"]
        prediction = predictors[target_id].predict(smiles)
        score = _target_score(prediction)
        weight = float(target["weight"])
        predictions.append({
            "target_chembl_id": target_id,
            "target_name": target.get("target_name") or target_id,
            "weight": weight,
            "direction": target.get("direction", "inhibit"),
            "target_score": score,
            **prediction,
        })
        weighted_score += weight * score
        weight_total += weight
    aggregate = weighted_score / weight_total
    aggregate -= min(len(alerts), 4) * 0.15
    all_in_domain = all(item["applicability_domain"]["inside"] for item in predictions)
    row = {
        "standard_inchi_key": compound["standard_inchi_key"],
        "canonical_smiles": smiles,
        "names": compound.get("names"),
        "reported_source_organisms": compound.get("organisms"),
        "natural_product_provenance": compound.get("provenance_ids"),
        "disease_id": specification["disease"].get("id") or specification["disease"].get("name"),
        "disease_name": specification["disease"].get("name") or specification["disease"].get("id"),
        "prioritization_score": round(aggregate, 6),
        "all_targets_inside_domain": all_in_domain,
        "structural_alerts": alerts,
        "descriptors": descriptors,
        "target_predictions": predictions,
        "research_status": "PRIORITIZED_FOR_EXPERIMENTAL_VALIDATION" if all_in_domain else "OUT_OF_DOMAIN_REVIEW_REQUIRED",
    }
    return row, None


def _write_results(output_dir: Path, rows: list[dict[str, Any]]) -> None:
    jsonl_path = output_dir / "ranked_candidates.jsonl"
    csv_path = output_dir / "ranked_candidates.csv"
    with jsonl_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
    columns = [
        "rank", "standard_inchi_key", "names", "canonical_smiles", "reported_source_organisms",
        "natural_product_provenance", "disease_id", "disease_name", "prioritization_score",
        "all_targets_inside_domain", "research_status", "structural_alerts", "target_predictions", "descriptors",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                **row,
                "structural_alerts": json.dumps(row["structural_alerts"], ensure_ascii=False),
                "target_predictions": json.dumps(row["target_predictions"], ensure_ascii=False),
                "descriptors": json.dumps(row["descriptors"], ensure_ascii=False),
            })


def _enrich_graph(graph_path: Path, rows: list[dict[str, Any]], specification: Mapping[str, Any]) -> None:
    disease_id = str(specification["disease"].get("id") or specification["disease"].get("name"))
    with KnowledgeGraph(graph_path) as graph:
        disease_node = graph.add_node("disease", disease_id, name=specification["disease"].get("name") or disease_id)
        for row in rows:
            compound_node = graph.add_node(
                "compound",
                row["standard_inchi_key"],
                name=row.get("names") or row["standard_inchi_key"],
                properties={"canonical_smiles": row["canonical_smiles"]},
            )
            graph.add_edge(
                compound_node,
                "PRIORITIZED_FOR_DISEASE",
                disease_node,
                evidence_source="PhytoTruth-screening",
                evidence_id=f"{disease_id}:rank:{row['rank']}",
                properties={"rank": row["rank"], "score": row["prioritization_score"], "status": row["research_status"]},
            )
            for prediction in row["target_predictions"]:
                target_id = prediction["target_chembl_id"]
                target_node = graph.add_node("target", target_id, name=prediction["target_name"])
                graph.add_edge(
                    compound_node,
                    "PREDICTED_AGAINST",
                    target_node,
                    evidence_source="PhytoTruth-ML",
                    evidence_id=f"{disease_id}:{target_id}:{row['standard_inchi_key']}",
                    properties={
                        "predicted_pchembl": prediction["predicted_pchembl"],
                        "target_score": prediction["target_score"],
                        "inside_applicability_domain": prediction["applicability_domain"]["inside"],
                        "prediction_interval": prediction["prediction_interval"],
                    },
                )
        graph.connection.commit()


def screen_natural_products(
    specifications: Path,
    natural_product_database: Path,
    output_dir: Path,
    *,
    knowledge_graph: Path | None = None,
    source: str | None = None,
    max_compounds: int | None = None,
    top_k: int = 100,
    out_of_domain_quota: int = 0,
    predictors: Mapping[str, Predictor] | None = None,
) -> dict[str, Any]:
    """Stream a natural-product index and retain the best auditable multi-target candidates."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    if out_of_domain_quota < 0 or out_of_domain_quota >= top_k:
        raise ValueError("out_of_domain_quota must be between 0 and top_k - 1")
    specifications = Path(specifications)
    natural_product_database = Path(natural_product_database)
    output_dir = Path(output_dir)
    if not natural_product_database.is_file():
        raise FileNotFoundError(natural_product_database)
    specification = _load_specification(specifications)
    predictor_map: Mapping[str, Predictor] = predictors or {
        target["target_chembl_id"]: load_predictor(str(target["model_dir"]))
        for target in specification["targets"]
    }
    missing = {target["target_chembl_id"] for target in specification["targets"]} - set(predictor_map)
    if missing:
        raise ValueError(f"Missing predictors for targets: {sorted(missing)}")

    heap: list[tuple[float, int, dict[str, Any]]] = []
    out_of_domain_heap: list[tuple[float, int, dict[str, Any]]] = []
    screened = property_filtered = prediction_failed = 0
    filter_reason_counts: dict[str, int] = {}
    prediction_error_samples: list[dict[str, str]] = []
    for sequence, compound in enumerate(
        _iter_compounds(natural_product_database, source=source, max_compounds=max_compounds), start=1
    ):
        screened += 1
        try:
            row, rejection = _candidate_row(compound, specification, predictor_map)
        except Exception as exc:
            prediction_failed += 1
            if len(prediction_error_samples) < 20:
                prediction_error_samples.append({
                    "standard_inchi_key": str(compound.get("standard_inchi_key")),
                    "error_type": type(exc).__name__,
                    "message": str(exc)[:500],
                })
            continue
        if rejection:
            property_filtered += 1
            for reason in rejection.split(","):
                filter_reason_counts[reason] = filter_reason_counts.get(reason, 0) + 1
            continue
        assert row is not None
        item = (float(row["prioritization_score"]), sequence, row)
        if len(heap) < top_k:
            heapq.heappush(heap, item)
        elif item[0] > heap[0][0]:
            heapq.heapreplace(heap, item)

        if out_of_domain_quota and not row["all_targets_inside_domain"]:
            if len(out_of_domain_heap) < out_of_domain_quota:
                heapq.heappush(out_of_domain_heap, item)
            elif item[0] > out_of_domain_heap[0][0]:
                heapq.heapreplace(out_of_domain_heap, item)

    rows = [item[2] for item in sorted(heap, key=lambda value: (-value[0], value[1]))]
    current_out_of_domain = sum(not row["all_targets_inside_domain"] for row in rows)
    if current_out_of_domain < out_of_domain_quota:
        selected_keys = {row["standard_inchi_key"] for row in rows}
        reserve = [item[2] for item in sorted(out_of_domain_heap, key=lambda value: (-value[0], value[1]))]
        for candidate in reserve:
            if current_out_of_domain >= out_of_domain_quota:
                break
            if candidate["standard_inchi_key"] in selected_keys:
                continue
            replace_index = next(
                (index for index in range(len(rows) - 1, -1, -1) if rows[index]["all_targets_inside_domain"]),
                None,
            )
            if replace_index is None:
                break
            removed = rows.pop(replace_index)
            selected_keys.remove(removed["standard_inchi_key"])
            rows.append(candidate)
            selected_keys.add(candidate["standard_inchi_key"])
            current_out_of_domain += 1
        rows.sort(key=lambda row: -float(row["prioritization_score"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_results(output_dir, rows)
    if knowledge_graph is not None:
        _enrich_graph(Path(knowledge_graph), rows, specification)
    summary = {
        "created_at": datetime.now(UTC).isoformat(),
        "screened_compounds": screened,
        "property_filtered_compounds": property_filtered,
        "property_filter_reason_counts": filter_reason_counts,
        "prediction_failed_compounds": prediction_failed,
        "prediction_error_samples": prediction_error_samples,
        "ranked_candidates": len(rows),
        "top_k": top_k,
        "out_of_domain_quota": out_of_domain_quota,
        "in_domain_candidates": sum(row["all_targets_inside_domain"] for row in rows),
        "out_of_domain_candidates": sum(not row["all_targets_inside_domain"] for row in rows),
        "source_filter": source,
        "max_compounds": max_compounds,
        "disease": specification["disease"],
        "targets": [
            {key: value for key, value in target.items() if key != "model_dir"} | {"model_dir": str(target["model_dir"])}
            for target in specification["targets"]
        ],
        "filters": specification["filters"],
        "score_semantics": (
            "weighted predicted pChEMBL adjusted for applicability-domain membership, nearest-training "
            "similarity, interval width, and structural-alert burden"
        ),
        "limitations": [
            "Predictions prioritize experimental follow-up; they do not establish efficacy, safety, or clinical utility.",
            "PAINS/Brenk matches are research alerts, not toxicity predictions.",
            "Reported source organisms are upstream provenance and are not independently taxonomically validated.",
        ],
        "specifications": str(specifications),
        "specifications_sha256": _sha256(specifications),
        "natural_product_database": str(natural_product_database),
        "natural_product_database_sha256": _sha256(natural_product_database),
        "knowledge_graph": str(knowledge_graph) if knowledge_graph else None,
        "outputs": {
            "jsonl": str(output_dir / "ranked_candidates.jsonl"),
            "csv": str(output_dir / "ranked_candidates.csv"),
        },
    }
    (output_dir / "screening_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary

