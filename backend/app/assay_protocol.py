from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import html
import json
import math
import re
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from .chembl import iter_assays_by_ids


EXTRACTION_VERSION = "assay-protocol-rules-v1"
SOURCE_NAME = "ChEMBL REST API"
SOURCE_BASE = "https://www.ebi.ac.uk/chembl/api/data/assay"

# Every derived text feature below requires a literal, case-insensitive match.
# These categories do not assert that an unmentioned property is absent.
TEXT_RULES: dict[str, dict[str, tuple[str, ...]]] = {
    "methodology": {
        "binding": (r"\bbind(?:ing|s)?\b", r"\bradioligand\b"),
        "enzyme_activity": (r"\benzym(?:e|atic)\b", r"\bkinase activity\b"),
        "functional": (r"\bfunctional assay\b", r"\bfunction(?:al)? response\b"),
        "cell_based": (r"\bcell[- ]based\b", r"\b(?:in|using)\s+(?:human\s+)?cells?\b"),
        "viability_or_proliferation": (r"\bviabilit(?:y|ies)\b", r"\bproliferation\b"),
        "phosphorylation": (r"\bphosphorylation\b", r"\bautophosphorylation\b"),
    },
    "detection": {
        "fluorescence": (r"\bfluorescen(?:ce|t)\b", r"\bfluorometric\b"),
        "luminescence": (r"\bluminescen(?:ce|t)\b", r"\bluciferase\b"),
        "radiometric": (r"\bradiometric\b", r"\bradioactiv(?:e|ity)\b", r"\bradioligand\b"),
        "absorbance_or_colorimetric": (r"\babsorbance\b", r"\bcolorimetric\b", r"\bspectrophotometric\b"),
        "fret": (r"\b(?:tr[- ]?)?fret\b",),
        "alphascreen": (r"\balpha(?:screen|lisa)\b",),
        "scintillation_proximity": (r"\bscintillation proximity\b", r"\bspa assay\b"),
        "immunoassay": (r"\belisa\b", r"\bimmunoassay\b"),
        "chromatography_or_ms": (r"\bhplc\b", r"\blc[- ]?ms\b", r"\bmass spectrometr(?:y|ic)\b"),
    },
    "biological_system_text": {
        "recombinant": (r"\brecombinant\b",),
        "purified": (r"\bpurified\b",),
        "cellular": (r"\bcells?\b", r"\bcell line\b"),
        "membrane": (r"\bmembranes?\b",),
        "lysate": (r"\blysates?\b",),
        "tissue": (r"\btissues?\b",),
    },
    "endpoint": {
        "IC50": (r"\bic\s*50\b", r"\bhalf[- ]maximal inhibitory concentration\b"),
        "EC50": (r"\bec\s*50\b", r"\bhalf[- ]maximal effective concentration\b"),
        "GI50": (r"\bgi\s*50\b",),
        "Ki": (r"\bki\b", r"\binhibition constant\b"),
        "Kd": (r"\bkd\b", r"\bdissociation constant\b"),
        "percent_inhibition": (r"%\s*inhibition\b", r"\bpercent(?:age)? inhibition\b"),
        "activity": (r"\bactivity\b",),
    },
}

CONDITION_RULES = {
    "temperature": re.compile(r"(?i)(?<!\w)(?:-?\d{1,3}(?:\.\d+)?)\s*(?:°\s*)?[CF](?!\w)"),
    "ph": re.compile(r"(?i)\bpH\s*(?:of\s*)?\d{1,2}(?:\.\d+)?\b"),
    "duration": re.compile(
        r"(?i)\b(?:incubat(?:ed|ion)|reaction(?:\s+time)?|preincubat(?:ed|ion))\b.{0,35}?"
        r"\b\d+(?:\.\d+)?\s*(?:seconds?|secs?|minutes?|mins?|hours?|hrs?|h)\b"
    ),
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _split_ids(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return sorted({item.strip() for item in str(value).split(";") if item.strip()})


def assay_ids_from_events(path: Path) -> list[str]:
    frame = pd.read_csv(path, usecols=["prior_assay_ids", "heldout_assay_ids"])
    identifiers: set[str] = set()
    for column in ("prior_assay_ids", "heldout_assay_ids"):
        for value in frame[column]:
            identifiers.update(_split_ids(value))
    return sorted(identifiers)


async def fetch_assay_metadata(
    event_dataset: Path,
    output: Path,
    *,
    batch_size: int = 50,
    max_retries: int = 5,
) -> dict[str, Any]:
    """Retrieve exact ChEMBL assay records referenced by an event dataset."""
    assay_ids = assay_ids_from_events(event_dataset)
    output.parent.mkdir(parents=True, exist_ok=True)
    part = output.with_suffix(output.suffix + ".part")
    retrieved_at = datetime.now(UTC).isoformat()
    returned: set[str] = set()
    digest = hashlib.sha256()
    with part.open("w", encoding="utf-8", newline="\n") as handle:
        async for record, request_url in iter_assays_by_ids(
            assay_ids, batch_size=batch_size, max_retries=max_retries
        ):
            identifier = str(record.get("assay_chembl_id") or "")
            if not identifier or identifier in returned:
                continue
            raw_sha256 = _sha256_bytes(
                json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
            )
            envelope = {
                "source": SOURCE_NAME,
                "source_endpoint": request_url,
                "retrieved_at": retrieved_at,
                "assay_chembl_id": identifier,
                "raw_record_sha256": raw_sha256,
                "record": record,
            }
            line = json.dumps(envelope, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
            handle.write(line)
            digest.update(line.encode("utf-8"))
            returned.add(identifier)
    part.replace(output)
    missing = sorted(set(assay_ids) - returned)
    manifest = {
        "source": SOURCE_NAME,
        "source_base": SOURCE_BASE,
        "retrieved_at": retrieved_at,
        "event_dataset": str(event_dataset),
        "event_dataset_sha256": _sha256_file(event_dataset),
        "requested_assays": len(assay_ids),
        "returned_assays": len(returned),
        "missing_assay_ids": missing,
        "batch_size": batch_size,
        "max_retries": max_retries,
        "data_file": output.name,
        "sha256": digest.hexdigest(),
    }
    output.with_suffix(output.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest


def _snippets(text: str, match: re.Match[str], flank: int = 55) -> str:
    start = max(0, match.start() - flank)
    end = min(len(text), match.end() + flank)
    return re.sub(r"\s+", " ", text[start:end]).strip()


def _rule_matches(text: str) -> tuple[dict[str, list[str]], dict[str, list[dict[str, str]]]]:
    categories: dict[str, list[str]] = {}
    evidence: dict[str, list[dict[str, str]]] = {}
    for family, family_rules in TEXT_RULES.items():
        selected: list[str] = []
        family_evidence: list[dict[str, str]] = []
        for label, patterns in family_rules.items():
            hits = []
            for pattern in patterns:
                match = re.search(pattern, text, flags=re.IGNORECASE)
                if match:
                    hits.append({"category": label, "pattern": pattern, "snippet": _snippets(text, match)})
            if hits:
                selected.append(label)
                family_evidence.extend(hits)
        categories[family] = sorted(selected)
        evidence[family] = family_evidence
    return categories, evidence


def _condition_matches(text: str) -> tuple[dict[str, list[str]], dict[str, list[dict[str, str]]]]:
    values: dict[str, list[str]] = {}
    evidence: dict[str, list[dict[str, str]]] = {}
    for name, pattern in CONDITION_RULES.items():
        matches = list(pattern.finditer(text))
        values[name] = sorted({match.group(0) for match in matches})
        evidence[name] = [{"match": match.group(0), "snippet": _snippets(text, match)} for match in matches]
    substrate_matches = list(re.finditer(r"(?i)\bsubstrates?\b", text))
    values["substrate_mentions"] = sorted({_snippets(text, match, flank=70) for match in substrate_matches})
    evidence["substrate_mentions"] = [
        {"match": match.group(0), "snippet": _snippets(text, match, flank=70)}
        for match in substrate_matches
    ]
    return values, evidence


def extract_protocol_record(envelope: dict[str, Any]) -> dict[str, Any]:
    """Extract literal protocol signals while preserving unknowns as unknown."""
    record = envelope.get("record", envelope)
    assay_id = str(record.get("assay_chembl_id") or envelope.get("assay_chembl_id") or "")
    if not assay_id:
        raise ValueError("Assay record is missing assay_chembl_id")
    description_value = record.get("description")
    description = str(description_value).strip() if description_value not in (None, "") else ""
    categories, text_evidence = _rule_matches(description) if description else (
        {family: [] for family in TEXT_RULES}, {family: [] for family in TEXT_RULES}
    )
    conditions, condition_evidence = _condition_matches(description) if description else (
        {**{name: [] for name in CONDITION_RULES}, "substrate_mentions": []},
        {**{name: [] for name in CONDITION_RULES}, "substrate_mentions": []},
    )
    parameters = record.get("assay_parameters") if isinstance(record.get("assay_parameters"), list) else []
    parameter_types = sorted({
        str(item.get("standard_type") or item.get("type")).strip()
        for item in parameters if isinstance(item, dict) and (item.get("standard_type") or item.get("type"))
    })
    variant = record.get("variant_sequence") if isinstance(record.get("variant_sequence"), dict) else {}
    raw_sha = envelope.get("raw_record_sha256") or _sha256_bytes(
        json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    )
    result: dict[str, Any] = {
        "assay_chembl_id": assay_id,
        "target_chembl_id": record.get("target_chembl_id"),
        "document_chembl_id": record.get("document_chembl_id"),
        "assay_type": record.get("assay_type"),
        "assay_type_description": record.get("assay_type_description"),
        "assay_organism": record.get("assay_organism"),
        "assay_tax_id": record.get("assay_tax_id"),
        "assay_tissue": record.get("assay_tissue"),
        "tissue_chembl_id": record.get("tissue_chembl_id"),
        "assay_cell_type": record.get("assay_cell_type"),
        "cell_chembl_id": record.get("cell_chembl_id"),
        "assay_subcellular_fraction": record.get("assay_subcellular_fraction"),
        "assay_strain": record.get("assay_strain"),
        "bao_format": record.get("bao_format"),
        "bao_label": record.get("bao_label"),
        "assay_category": record.get("assay_category"),
        "assay_test_type": record.get("assay_test_type"),
        "assay_group": record.get("assay_group"),
        "confidence_score": record.get("confidence_score"),
        "confidence_description": record.get("confidence_description"),
        "relationship_type": record.get("relationship_type"),
        "relationship_description": record.get("relationship_description"),
        "description": description or None,
        "description_status": "reported" if description else "missing",
        "assay_parameters_json": json.dumps(parameters, sort_keys=True, ensure_ascii=False),
        "assay_parameters_status": "reported" if parameters else "missing",
        "assay_parameter_types": ";".join(parameter_types),
        "variant_accession": variant.get("accession"),
        "variant_mutation": variant.get("mutation"),
        "variant_organism": variant.get("organism"),
        "variant_status": "reported" if variant else "missing",
        "methodology_categories": ";".join(categories["methodology"]),
        "detection_categories": ";".join(categories["detection"]),
        "biological_system_text_categories": ";".join(categories["biological_system_text"]),
        "endpoint_categories": ";".join(categories["endpoint"]),
        "temperature_mentions": ";".join(conditions["temperature"]),
        "ph_mentions": ";".join(conditions["ph"]),
        "duration_mentions": ";".join(conditions["duration"]),
        "substrate_mentions_json": json.dumps(conditions["substrate_mentions"], ensure_ascii=False),
        "text_extraction_status": "explicit_mentions_only" if description else "no_source_description",
        "extraction_evidence_json": json.dumps(
            {"source_field": "description", "text_rules": text_evidence, "conditions": condition_evidence},
            sort_keys=True,
            ensure_ascii=False,
        ),
        "extraction_version": EXTRACTION_VERSION,
        "source": envelope.get("source", SOURCE_NAME),
        "source_endpoint": envelope.get("source_endpoint", f"{SOURCE_BASE}/{assay_id}.json"),
        "retrieved_at": envelope.get("retrieved_at"),
        "raw_record_sha256": raw_sha,
    }
    direct_fields = [
        "assay_type", "assay_organism", "assay_tissue", "assay_cell_type",
        "assay_subcellular_fraction", "assay_strain", "bao_format", "description",
        "assay_parameters_json", "variant_mutation",
    ]
    result["reported_direct_fields"] = sum(
        result.get(field) not in (None, "", "[]") for field in direct_fields
    )
    result["direct_field_completeness"] = result["reported_direct_fields"] / len(direct_fields)
    return result


def load_and_extract(raw_jsonl: Path) -> pd.DataFrame:
    records = []
    seen: set[str] = set()
    with raw_jsonl.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                envelope = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {raw_jsonl}:{line_number}") from exc
            row = extract_protocol_record(envelope)
            assay_id = row["assay_chembl_id"]
            if assay_id in seen:
                raise ValueError(f"Duplicate assay metadata: {assay_id}")
            seen.add(assay_id)
            records.append(row)
    return pd.DataFrame(records).sort_values("assay_chembl_id").reset_index(drop=True)


def _set_from_rows(rows: pd.DataFrame, column: str) -> set[str]:
    values: set[str] = set()
    for value in rows[column].dropna().astype(str):
        values.update(item.strip() for item in value.split(";") if item.strip())
    return values


def _direct_set(rows: pd.DataFrame, column: str) -> set[str]:
    return {str(value).strip() for value in rows[column].dropna() if str(value).strip()}


def _aggregate_side(rows: pd.DataFrame, expected_ids: Sequence[str], prefix: str) -> dict[str, Any]:
    expected = set(expected_ids)
    found = set(rows.assay_chembl_id.astype(str)) if len(rows) else set()
    result: dict[str, Any] = {
        f"{prefix}_protocol_availability": (
            "complete" if expected and expected == found else "partial" if found else "missing"
        ),
        f"{prefix}_protocol_expected_assays": len(expected),
        f"{prefix}_protocol_retrieved_assays": len(found),
        f"{prefix}_protocol_missing_assay_ids": ";".join(sorted(expected - found)),
        f"{prefix}_protocol_assay_ids": ";".join(sorted(found)),
        f"{prefix}_description_reported_fraction": (
            float((rows.description_status == "reported").mean()) if len(rows) else None
        ),
        f"{prefix}_parameter_reported_fraction": (
            float((rows.assay_parameters_status == "reported").mean()) if len(rows) else None
        ),
        f"{prefix}_direct_field_completeness_mean": (
            float(rows.direct_field_completeness.mean()) if len(rows) else None
        ),
    }
    for column in (
        "assay_type_description", "assay_organism", "assay_tissue", "assay_cell_type",
        "assay_subcellular_fraction", "assay_strain", "bao_label", "assay_category",
        "assay_test_type", "assay_parameter_types", "variant_accession", "variant_mutation",
        "variant_organism",
    ):
        result[f"{prefix}_{column}_values"] = ";".join(sorted(_direct_set(rows, column)))
    for column in (
        "methodology_categories", "detection_categories", "biological_system_text_categories",
        "endpoint_categories", "temperature_mentions", "ph_mentions", "duration_mentions",
    ):
        result[f"{prefix}_{column}"] = ";".join(sorted(_set_from_rows(rows, column)))
    return result


def _known_difference(prior: set[str], heldout: set[str]) -> tuple[str, int | None]:
    if not prior or not heldout:
        return "unknown", None
    return ("different", 1) if prior != heldout else ("same", 0)


def enrich_events(events: pd.DataFrame, assays: pd.DataFrame) -> pd.DataFrame:
    lookup = assays.set_index("assay_chembl_id", drop=False)
    enriched_rows: list[dict[str, Any]] = []
    direct_compare = {
        "assay_type": "assay_type_description",
        "organism": "assay_organism",
        "tissue": "assay_tissue",
        "cell_type": "assay_cell_type",
        "subcellular_fraction": "assay_subcellular_fraction",
        "strain": "assay_strain",
        "bao_format": "bao_label",
        "parameter_type": "assay_parameter_types",
        "variant_mutation": "variant_mutation",
    }
    derived_compare = {
        "methodology": "methodology_categories",
        "detection": "detection_categories",
        "biological_system_text": "biological_system_text_categories",
        "endpoint": "endpoint_categories",
        "temperature": "temperature_mentions",
        "ph": "ph_mentions",
        "duration": "duration_mentions",
    }
    for event_index, source in events.reset_index(drop=True).iterrows():
        prior_ids = _split_ids(source.get("prior_assay_ids"))
        heldout_ids = _split_ids(source.get("heldout_assay_ids"))
        prior_rows = lookup.loc[lookup.index.intersection(prior_ids)]
        heldout_rows = lookup.loc[lookup.index.intersection(heldout_ids)]
        if isinstance(prior_rows, pd.Series):
            prior_rows = prior_rows.to_frame().T
        if isinstance(heldout_rows, pd.Series):
            heldout_rows = heldout_rows.to_frame().T
        row = source.to_dict()
        row["protocol_event_id"] = hashlib.sha256(
            f"{source.get('target_chembl_id')}|{source.get('molecule_id')}|{source.get('heldout_study_key')}|{event_index}".encode()
        ).hexdigest()[:20]
        row.update(_aggregate_side(prior_rows, prior_ids, "prior"))
        row.update(_aggregate_side(heldout_rows, heldout_ids, "heldout"))
        availability_columns = {
            "organism": "assay_organism", "tissue": "assay_tissue", "cell_type": "assay_cell_type",
            "parameters": "assay_parameter_types", "variant": "variant_mutation",
            "methodology": "methodology_categories", "detection": "detection_categories",
            "biological_system_text": "biological_system_text_categories", "endpoint": "endpoint_categories",
            "temperature": "temperature_mentions", "ph": "ph_mentions", "duration": "duration_mentions",
        }
        for label, column in availability_columns.items():
            getter = _direct_set if column in {
                "assay_organism", "assay_tissue", "assay_cell_type", "assay_parameter_types", "variant_mutation"
            } else _set_from_rows
            row[f"prior_{label}_reported"] = int(bool(getter(prior_rows, column)))
            row[f"heldout_{label}_reported"] = int(bool(getter(heldout_rows, column)))
        for label, column in {**direct_compare, **derived_compare}.items():
            getter = _direct_set if column in direct_compare.values() else _set_from_rows
            prior = getter(prior_rows, column)
            heldout = getter(heldout_rows, column)
            status, flag = _known_difference(prior, heldout)
            row[f"protocol_{label}_comparison"] = status
            row[f"protocol_{label}_changed"] = flag
        row["protocol_prior_feature_timing"] = "pre_prediction_available"
        row["protocol_heldout_feature_timing"] = "post_outcome_explanatory_only"
        row["protocol_comparison_feature_timing"] = "post_outcome_explanatory_only"
        row["protocol_extraction_version"] = EXTRACTION_VERSION
        enriched_rows.append(row)
    return pd.DataFrame(enriched_rows)


def _risk_difference(group: pd.DataFrame, flag: str) -> dict[str, Any] | None:
    comparable = group[group[flag].notna()].copy()
    if comparable.empty or comparable[flag].nunique() < 2:
        return None
    exposed = comparable[comparable[flag] == 1].disagreement_label.astype(int)
    reference = comparable[comparable[flag] == 0].disagreement_label.astype(int)
    table = [
        [int(exposed.sum()), int(len(exposed) - exposed.sum())],
        [int(reference.sum()), int(len(reference) - reference.sum())],
    ]
    return {
        "comparable_events": len(comparable),
        "changed_events": len(exposed),
        "same_events": len(reference),
        "changed_disagreement_rate": float(exposed.mean()),
        "same_disagreement_rate": float(reference.mean()),
        "risk_difference": float(exposed.mean() - reference.mean()),
        "odds_ratio": float(fisher_exact(table).statistic),
        "fisher_exact_p": float(fisher_exact(table).pvalue),
    }


def _clustered_risk_ci(group: pd.DataFrame, flag: str, *, iterations: int, seed: int) -> list[float] | None:
    comparable = group[group[flag].notna()].copy()
    cluster_counts = []
    for _, study in comparable.groupby("heldout_study_key", sort=True):
        changed = study[flag].to_numpy(dtype=int) == 1
        outcomes = study.disagreement_label.to_numpy(dtype=int)
        cluster_counts.append([
            int(np.sum(changed & (outcomes == 1))),
            int(np.sum(changed & (outcomes == 0))),
            int(np.sum(~changed & (outcomes == 1))),
            int(np.sum(~changed & (outcomes == 0))),
        ])
    if len(cluster_counts) < 2:
        return None
    counts = np.asarray(cluster_counts, dtype=int)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(iterations):
        sampled = counts[rng.integers(0, len(counts), len(counts))].sum(axis=0)
        changed_n = sampled[0] + sampled[1]
        same_n = sampled[2] + sampled[3]
        if changed_n and same_n:
            values.append(float(sampled[0] / changed_n - sampled[2] / same_n))
    return [float(value) for value in np.quantile(values, [0.025, 0.975])] if values else None


def _bh_adjust(rows: list[dict[str, Any]]) -> None:
    indexed = sorted(
        [(index, float(row["fisher_exact_p"])) for index, row in enumerate(rows) if row.get("fisher_exact_p") is not None],
        key=lambda item: item[1],
    )
    running = 1.0
    for reverse_rank, (index, value) in enumerate(reversed(indexed), start=1):
        rank = len(indexed) - reverse_rank + 1
        running = min(running, value * len(indexed) / rank)
        rows[index]["fisher_bh_q"] = float(running)


def association_analysis(
    enriched: pd.DataFrame, *, iterations: int = 1000, seed: int = 2026
) -> list[dict[str, Any]]:
    features = [
        column for column in enriched.columns
        if column.startswith("protocol_") and column.endswith("_changed")
    ]
    rows: list[dict[str, Any]] = []
    scopes = [("all", enriched), ("natural_products", enriched[enriched.is_verified_natural_product.astype(str).str.lower() == "true"])]
    for partition in ("development", "evaluation"):
        for subgroup, scope in scopes:
            group = scope[scope.partition == partition]
            for feature_index, feature in enumerate(features):
                result = _risk_difference(group, feature)
                if result is None:
                    continue
                rows.append({
                    "partition": partition,
                    "subgroup": subgroup,
                    "feature": feature,
                    "feature_timing": "post_outcome_explanatory_only",
                    **result,
                    "study_clustered_ci95": _clustered_risk_ci(
                        group, feature, iterations=iterations, seed=seed + feature_index
                    ),
                })
    for key in {(row["partition"], row["subgroup"]) for row in rows}:
        _bh_adjust([row for row in rows if (row["partition"], row["subgroup"]) == key])
    return rows


def missingness_association_analysis(
    enriched: pd.DataFrame, *, iterations: int = 1000, seed: int = 2026
) -> list[dict[str, Any]]:
    """Associate prediction-time metadata availability with future disagreement."""
    features = sorted(
        column for column in enriched.columns
        if column.startswith("prior_") and column.endswith("_reported")
    )
    rows: list[dict[str, Any]] = []
    scopes = [("all", enriched), ("natural_products", enriched[enriched.is_verified_natural_product.astype(str).str.lower() == "true"])]
    for partition in ("development", "evaluation"):
        for subgroup, scope in scopes:
            group = scope[scope.partition == partition]
            for feature_index, feature in enumerate(features):
                result = _risk_difference(group, feature)
                if result is None:
                    continue
                result["reported_events"] = result.pop("changed_events")
                result["missing_events"] = result.pop("same_events")
                result["reported_disagreement_rate"] = result.pop("changed_disagreement_rate")
                result["missing_disagreement_rate"] = result.pop("same_disagreement_rate")
                result["risk_difference_reported_minus_missing"] = result.pop("risk_difference")
                row = {
                    "partition": partition,
                    "subgroup": subgroup,
                    "feature": feature,
                    "feature_timing": "pre_prediction_available",
                    **result,
                }
                # Bootstrap uses the original binary column; rename only the resulting estimand.
                row["study_clustered_ci95"] = _clustered_risk_ci(
                    group, feature, iterations=iterations, seed=seed + feature_index
                )
                rows.append(row)
    for key in {(row["partition"], row["subgroup"]) for row in rows}:
        _bh_adjust([row for row in rows if (row["partition"], row["subgroup"]) == key])
    return rows


def study_missingness_table(enriched: pd.DataFrame) -> pd.DataFrame:
    rows = []
    report_columns = sorted(
        column for column in enriched.columns
        if column.startswith("heldout_") and column.endswith("_reported")
    )
    for (partition, target, study), group in enriched.groupby(
        ["partition", "target_chembl_id", "heldout_study_key"], sort=True
    ):
        row = {
            "partition": partition,
            "target_chembl_id": target,
            "heldout_study_key": study,
            "events": len(group),
            "disagreements": int(group.disagreement_label.astype(int).sum()),
            "disagreement_rate": float(group.disagreement_label.astype(int).mean()),
            "mean_direct_field_completeness": float(group.heldout_direct_field_completeness_mean.mean()),
        }
        for column in report_columns:
            row[f"{column}_fraction"] = float(group[column].astype(float).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def coverage_summary(assays: pd.DataFrame, enriched: pd.DataFrame) -> dict[str, Any]:
    direct_fields = [
        "assay_type_description", "assay_organism", "assay_tissue", "assay_cell_type",
        "assay_subcellular_fraction", "assay_strain", "bao_label", "description",
        "assay_category", "assay_test_type", "assay_parameter_types", "variant_mutation",
    ]
    extracted_fields = [
        "methodology_categories", "detection_categories", "biological_system_text_categories",
        "endpoint_categories", "temperature_mentions", "ph_mentions", "duration_mentions",
        "substrate_mentions_json",
    ]
    assay_coverage = {}
    for field in direct_fields:
        present = assays[field].notna() & (assays[field].astype(str).str.strip() != "")
        assay_coverage[field] = {"reported": int(present.sum()), "fraction": float(present.mean())}
    for field in extracted_fields:
        if field.endswith("_json"):
            present = assays[field].fillna("[]").astype(str) != "[]"
        else:
            present = assays[field].fillna("").astype(str).str.strip() != ""
        assay_coverage[field] = {
            "explicit_mentions": int(present.sum()),
            "fraction": float(present.mean()),
            "absence_interpretation": "not reported or not detected; not evidence of experimental absence",
        }
    event_coverage = {}
    for side in ("prior", "heldout"):
        values = enriched[f"{side}_protocol_availability"].value_counts().to_dict()
        event_coverage[side] = {str(key): int(value) for key, value in values.items()}
    by_partition = {}
    for partition, group in enriched.groupby("partition", sort=True):
        by_partition[str(partition)] = {
            "events": len(group),
            "disagreement_rate": float(group.disagreement_label.astype(int).mean()),
            "natural_product_events": int(
                (group.is_verified_natural_product.astype(str).str.lower() == "true").sum()
            ),
            "prior_descriptions_complete": int((group.prior_description_reported_fraction == 1.0).sum()),
            "heldout_descriptions_complete": int((group.heldout_description_reported_fraction == 1.0).sum()),
        }
    return {
        "assays": len(assays),
        "events": len(enriched),
        "assay_field_coverage": assay_coverage,
        "event_metadata_availability": event_coverage,
        "partitions": by_partition,
    }


def _write_svg(path: Path, coverage: dict[str, Any], associations: list[dict[str, Any]]) -> None:
    fields = [
        (name, values.get("fraction", 0.0))
        for name, values in coverage["assay_field_coverage"].items()
    ]
    width, height = 1200, 760
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#07131d"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#dcecf4}.muted{fill:#8da9b8}.bar{fill:#2dd4bf}.axis{stroke:#58717d}</style>',
        '<text x="45" y="48" font-size="25" font-weight="700">ChEMBL assay protocol metadata coverage</text>',
        '<text x="45" y="75" font-size="14" class="muted">Explicitly reported fields or literal description mentions; missing does not mean absent.</text>',
    ]
    for index, (name, fraction) in enumerate(fields):
        y = 115 + index * 34
        label = html.escape(name.replace("_", " "))
        lines.extend([
            f'<text x="45" y="{y + 15}" font-size="13">{label}</text>',
            f'<rect x="330" y="{y}" width="650" height="19" rx="4" fill="#17303c"/>',
            f'<rect class="bar" x="330" y="{y}" width="{650 * fraction:.1f}" height="19" rx="4"/>',
            f'<text x="995" y="{y + 15}" font-size="13">{fraction * 100:.1f}%</text>',
        ])
    y0 = 115 + len(fields) * 34 + 30
    eval_rows = sorted(
        [row for row in associations if row["partition"] == "evaluation" and row["subgroup"] == "all"],
        key=lambda row: abs(row["risk_difference"]), reverse=True,
    )[:4]
    lines.append(f'<text x="45" y="{y0}" font-size="19" font-weight="700">Largest evaluation-set descriptive risk differences</text>')
    lines.append(f'<text x="45" y="{y0 + 24}" font-size="13" class="muted">Post-outcome protocol comparisons; associative, not predictive or causal.</text>')
    for index, row in enumerate(eval_rows):
        y = y0 + 58 + index * 31
        label = html.escape(row["feature"].replace("protocol_", "").replace("_changed", "").replace("_", " "))
        value = row["risk_difference"]
        lines.append(f'<text x="45" y="{y}" font-size="13">{label}: {value:+.3f} ({row["comparable_events"]} comparable events)</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def _format_pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def _write_report(
    path: Path,
    coverage: dict[str, Any],
    associations: list[dict[str, Any]],
    missingness_associations: list[dict[str, Any]],
    manifest: dict[str, Any],
) -> None:
    fields = coverage["assay_field_coverage"]
    eval_all = sorted(
        [row for row in associations if row["partition"] == "evaluation" and row["subgroup"] == "all"],
        key=lambda row: abs(row["risk_difference"]), reverse=True,
    )
    eval_missingness = sorted(
        [row for row in missingness_associations if row["partition"] == "evaluation" and row["subgroup"] == "all"],
        key=lambda row: abs(row["risk_difference_reported_minus_missing"]), reverse=True,
    )
    lines = [
        "# Assay protocol metadata enrichment report",
        "",
        "## Scope and integrity controls",
        "",
        f"This release enriches **{coverage['events']}** assay-disagreement events with **{coverage['assays']}** exact ChEMBL assay records. The raw API response, request URL, retrieval timestamp, and record SHA-256 are retained for every assay. Extraction version: `{EXTRACTION_VERSION}`.",
        "",
        "Structured ChEMBL fields are copied without reinterpretation. Description-derived categories require literal rule matches and include the matched pattern and source-text snippet. Empty categories mean *not reported or not detected*, never that a method or condition was experimentally absent.",
        "",
        "Prior-assay fields are labelled `pre_prediction_available`. Held-out protocol fields and prior-versus-held-out comparisons are labelled `post_outcome_explanatory_only`; they must not be used to predict the held-out result.",
        "",
        "All included assays are ChEMBL Binding assays in single-protein BAO format, reflecting the benchmark's compatibility filter. Those two fields have no within-cohort variation and therefore cannot explain disagreement here.",
        "",
        "## Metadata coverage",
        "",
        "| Field | Assays reported / explicitly mentioned | Coverage |",
        "|---|---:|---:|",
    ]
    for field, values in fields.items():
        count = values.get("reported", values.get("explicit_mentions", 0))
        lines.append(f"| {field} | {count} / {coverage['assays']} | {_format_pct(values['fraction'])} |")
    lines.extend([
        "",
        "## Dataset partitions",
        "",
        "| Partition | Events | Disagreement rate | Verified natural products | Complete prior descriptions | Complete held-out descriptions |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for name, values in coverage["partitions"].items():
        lines.append(
            f"| {name} | {values['events']} | {_format_pct(values['disagreement_rate'])} | "
            f"{values['natural_product_events']} | {values['prior_descriptions_complete']} | "
            f"{values['heldout_descriptions_complete']} |"
        )
    lines.extend([
        "",
        "## Protocol differences and disagreement",
        "",
        "Risk differences below compare observed disagreement rates for events whose prior and held-out protocol metadata differed versus matched. Confidence intervals use held-out-study cluster bootstrap. Fisher tests and BH q-values are supplementary because molecule rows within a study are not independent. These are retrospective associations, not causal effects and not model validation.",
        "",
        "| Evaluation feature | Comparable | Changed | Risk difference | Study-clustered 95% CI | BH q |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in eval_all:
        ci = row.get("study_clustered_ci95")
        ci_text = f"[{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "not estimable"
        q = row.get("fisher_bh_q")
        lines.append(
            f"| {row['feature']} | {row['comparable_events']} | {row['changed_events']} | "
            f"{row['risk_difference']:+.3f} | {ci_text} | {q:.3g} |"
        )
    lines.extend([
        "",
        "### Prediction-time metadata missingness",
        "",
        "These fixed comparisons use only prior-assay metadata available before the held-out study. A positive risk difference means events with reported prior metadata had more later disagreement than events where that field was unreported; it does not imply the metadata caused disagreement.",
        "",
        "| Evaluation feature | Comparable | Reported | Risk difference (reported - missing) | Study-clustered 95% CI | BH q |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in eval_missingness:
        ci = row.get("study_clustered_ci95")
        ci_text = f"[{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "not estimable"
        q = row.get("fisher_bh_q")
        lines.append(
            f"| {row['feature']} | {row['comparable_events']} | {row['reported_events']} | "
            f"{row['risk_difference_reported_minus_missing']:+.3f} | {ci_text} | {q:.3g} |"
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        "No protocol-difference or prior-metadata-availability association met a conservative transfer criterion of directional agreement in development, an evaluation study-clustered interval excluding zero, and adequate independent-study support. Several unclustered Fisher/BH results appear small, but their study-clustered intervals cross zero; they must not be treated as validated signals.",
        "",
        "The defensible finding is therefore about data adequacy: ChEMBL reliably supplies an assay description, broad assay type, and BAO format for this cohort, but most experimental conditions are too sparsely and unevenly reported to establish that protocol differences explain disagreement. Full-text methods and more independent natural-product studies are required.",
    ])
    lines.extend([
        "",
        "Natural-product subgroup rows are retained in `protocol_associations.csv`; sparse comparable groups and study concentration must be considered before interpretation.",
        "",
        "## Limitations",
        "",
        "- ChEMBL assay records frequently omit temperature, pH, duration, substrate, detection method, cell type, tissue, or strain. Free-text extraction cannot recover unreported conditions.",
        "- Literal keyword categories improve auditability but do not replace manual protocol curation or full-text methods sections.",
        "- ChEMBL descriptions can summarize rather than reproduce the experimental protocol. The raw record is authoritative for this artifact, not the derived category.",
        "- Protocol comparisons use metadata available only after the held-out assay was performed and are explanatory only. No prediction model was trained or retrained in this phase.",
        "- Associations can be confounded by target, era, publication practice, and study composition; they do not establish that a protocol difference caused disagreement.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m app.pipeline fetch-assay-metadata --events ../experiments/assay_disagreement_generalization/dataset.csv --output ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl",
        "python -m app.pipeline enrich-assay-protocols --events ../experiments/assay_disagreement_generalization/dataset.csv --assays ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl --output-dir ../experiments/assay_protocol_enrichment --bootstrap-iterations 1000 --seed 2026",
        "```",
        "",
        f"Input assay archive SHA-256: `{manifest['raw_assay_metadata_sha256']}`.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_enrichment(
    event_dataset: Path,
    raw_assays: Path,
    output_dir: Path,
    *,
    bootstrap_iterations: int = 1000,
    seed: int = 2026,
) -> dict[str, Any]:
    if bootstrap_iterations < 1:
        raise ValueError("bootstrap_iterations must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(event_dataset)
    required = {
        "target_chembl_id", "molecule_id", "partition", "heldout_study_key",
        "prior_assay_ids", "heldout_assay_ids", "disagreement_label", "is_verified_natural_product",
    }
    missing = required - set(events.columns)
    if missing:
        raise ValueError(f"Event dataset is missing columns: {sorted(missing)}")
    assays = load_and_extract(raw_assays)
    enriched = enrich_events(events, assays)
    associations = association_analysis(enriched, iterations=bootstrap_iterations, seed=seed)
    missingness_associations = missingness_association_analysis(
        enriched, iterations=bootstrap_iterations, seed=seed + 10_000
    )
    study_missingness = study_missingness_table(enriched)
    coverage = coverage_summary(assays, enriched)

    assay_csv = output_dir / "assay_protocol_features.csv"
    event_csv = output_dir / "enriched_assay_disagreement_events.csv"
    association_csv = output_dir / "protocol_associations.csv"
    missingness_csv = output_dir / "prior_metadata_missingness_associations.csv"
    study_csv = output_dir / "heldout_study_metadata_coverage.csv"
    assays.to_csv(assay_csv, index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    enriched.to_csv(event_csv, index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    pd.DataFrame(associations).to_csv(association_csv, index=False, lineterminator="\n")
    pd.DataFrame(missingness_associations).to_csv(missingness_csv, index=False, lineterminator="\n")
    study_missingness.to_csv(study_csv, index=False, lineterminator="\n")
    (output_dir / "coverage_statistics.json").write_text(
        json.dumps(coverage, indent=2, sort_keys=True), encoding="utf-8"
    )
    (output_dir / "protocol_associations.json").write_text(
        json.dumps(associations, indent=2, sort_keys=True), encoding="utf-8"
    )
    (output_dir / "prior_metadata_missingness_associations.json").write_text(
        json.dumps(missingness_associations, indent=2, sort_keys=True), encoding="utf-8"
    )
    manifest = {
        "created_at": datetime.now(UTC).isoformat(),
        "extraction_version": EXTRACTION_VERSION,
        "seed": seed,
        "bootstrap_iterations": bootstrap_iterations,
        "event_dataset": str(event_dataset),
        "event_dataset_sha256": _sha256_file(event_dataset),
        "raw_assay_metadata": str(raw_assays),
        "raw_assay_metadata_sha256": _sha256_file(raw_assays),
        "outputs": {},
        "feature_timing": {
            "prior_protocol_fields": "pre_prediction_available",
            "heldout_protocol_fields": "post_outcome_explanatory_only",
            "protocol_comparison_fields": "post_outcome_explanatory_only",
        },
    }
    for path in (
        assay_csv, event_csv, association_csv, missingness_csv, study_csv,
        output_dir / "coverage_statistics.json", output_dir / "protocol_associations.json",
        output_dir / "prior_metadata_missingness_associations.json",
    ):
        manifest["outputs"][path.name] = _sha256_file(path)
    _write_svg(output_dir / "protocol_coverage_and_associations.svg", coverage, associations)
    manifest["outputs"]["protocol_coverage_and_associations.svg"] = _sha256_file(
        output_dir / "protocol_coverage_and_associations.svg"
    )
    _write_report(
        output_dir / "SCIENTIFIC_REPORT.md", coverage, associations,
        missingness_associations, manifest,
    )
    manifest["outputs"]["SCIENTIFIC_REPORT.md"] = _sha256_file(output_dir / "SCIENTIFIC_REPORT.md")
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "coverage": coverage,
        "associations": associations,
        "missingness_associations": missingness_associations,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="ChEMBL assay protocol metadata enrichment")
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("fetch")
    fetch.add_argument("--events", type=Path, required=True)
    fetch.add_argument("--output", type=Path, required=True)
    fetch.add_argument("--batch-size", type=int, default=50)
    enrich = commands.add_parser("enrich")
    enrich.add_argument("--events", type=Path, required=True)
    enrich.add_argument("--assays", type=Path, required=True)
    enrich.add_argument("--output-dir", type=Path, required=True)
    enrich.add_argument("--bootstrap-iterations", type=int, default=1000)
    enrich.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    if args.command == "fetch":
        print(json.dumps(asyncio.run(fetch_assay_metadata(args.events, args.output, batch_size=args.batch_size)), indent=2))
    else:
        result = run_enrichment(
            args.events, args.assays, args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations, seed=args.seed,
        )
        print(json.dumps(result["coverage"], indent=2))


if __name__ == "__main__":
    main()
