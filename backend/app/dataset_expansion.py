from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .assay_disagreement import _jsonl, _sha256, build_dataset
from .assay_protocol import EXTRACTION_VERSION, enrich_events, load_and_extract
from .independent_validation import (
    assert_independence,
    build_previous_exposure_ledger,
    strict_independent_filter,
    validation_measurements,
)
from .multi_target_benchmark import load_target_specs


EXPANSION_VERSION = "assay-disagreement-dataset-expansion-v1"
DISAGREEMENT_THRESHOLD = 1.0


def merge_ledgers(*ledgers: dict[str, set[str]]) -> dict[str, set[str]]:
    keys = ("targets", "molecules", "activities", "assays", "documents")
    return {key: set().union(*(ledger[key] for ledger in ledgers)) for key in keys}


def resolve_target_specs(specifications: Path, data_root: Path, output: Path) -> Path:
    """Resolve frozen target identities to a local snapshot without altering selection."""
    resolved = []
    for spec in load_target_specs(specifications):
        target = str(spec["target_chembl_id"])
        item = dict(spec)
        item["activities"] = str((data_root / f"{target}.normalized.jsonl").resolve())
        item["natural_products"] = str((data_root / f"{target}.natural.jsonl").resolve())
        for key in ("activities", "natural_products"):
            if not Path(item[key]).is_file():
                raise FileNotFoundError(f"Missing frozen source snapshot: {item[key]}")
        resolved.append(item)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(resolved, indent=2), encoding="utf-8")
    return output


def _event_identifier(row: pd.Series | dict[str, Any]) -> str:
    value = row if isinstance(row, dict) else row.to_dict()
    payload = "|".join((
        str(value["target_chembl_id"]), str(value["molecule_id"]),
        str(value["heldout_study_key"]), str(value["heldout_year"]),
    ))
    return "ADE-" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")


def _write_jsonl(rows: Iterable[dict[str, Any]], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n")
            count += 1
    return count


def complete_measurement_records(events: pd.DataFrame, specs: Path) -> Iterable[dict[str, Any]]:
    roles: dict[str, tuple[str, str]] = {}
    for _, event in events.iterrows():
        event_id = str(event["expansion_event_id"])
        for role in ("prior", "heldout"):
            for activity in str(event[f"{role}_activity_ids"]).split(";"):
                if activity:
                    roles[activity] = (event_id, role)
    found: set[str] = set()
    for spec in load_target_specs(specs):
        for record in _jsonl(Path(spec["activities"])):
            activity = str(record.get("activity_id", record.get("source_record_id")))
            if activity in roles:
                found.add(activity)
                event_id, role = roles[activity]
                yield {
                    "expansion_event_id": event_id,
                    "evidence_role": role,
                    "target_chembl_id": spec["target_chembl_id"],
                    "normalized_record": record,
                }
    missing = set(roles) - found
    if missing:
        raise AssertionError(f"Used source activities are missing: {sorted(missing)[:5]}")


def natural_product_provenance(events: pd.DataFrame, specs: Path) -> pd.DataFrame:
    accepted = set(events.loc[events.is_verified_natural_product.astype(bool), "molecule_id"].astype(str))
    rows: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for spec in load_target_specs(specs):
        for record in _jsonl(Path(spec["natural_products"])):
            molecule = str(record.get("molecule_chembl_id"))
            if molecule not in accepted:
                continue
            structure = record.get("standardized_structure") or {}
            for source in record.get("natural_product_provenance") or []:
                key = (
                    molecule, str(source.get("source")), str(source.get("source_id")),
                    str(source.get("dataset_version")),
                )
                rows[key] = {
                    "molecule_id": molecule,
                    "target_chembl_id": spec["target_chembl_id"],
                    "source": source.get("source"),
                    "source_id": source.get("source_id"),
                    "dataset_version": source.get("dataset_version"),
                    "source_file": source.get("source_file"),
                    "source_row": source.get("source_row"),
                    "standard_inchi_key": structure.get("standard_inchi_key"),
                    "standard_inchi": structure.get("standard_inchi"),
                    "standardized_canonical_smiles": structure.get("canonical_smiles"),
                }
    columns = [
        "molecule_id", "target_chembl_id", "source", "source_id", "dataset_version",
        "source_file", "source_row", "standard_inchi_key", "standard_inchi",
        "standardized_canonical_smiles",
    ]
    return pd.DataFrame(sorted(rows.values(), key=lambda row: (row["molecule_id"], str(row["source"]), str(row["source_id"]))), columns=columns)


def prepare_expansion(
    specifications: Path,
    data_root: Path,
    previous_specs: list[Path],
    previous_events: list[Path],
    output_dir: Path,
    *,
    minimum_heldout_year: int = 2020,
) -> dict[str, Any]:
    if len(previous_specs) != len(previous_events):
        raise ValueError("Each previous specification requires its corresponding event dataset")
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = resolve_target_specs(specifications, data_root, output_dir / "resolved_targets.json")
    candidates, construction = build_dataset(resolved)
    ledgers = []
    exposure_summaries = []
    for specs, events in zip(previous_specs, previous_events, strict=True):
        ledger, summary = build_previous_exposure_ledger(specs, events)
        ledgers.append(ledger)
        exposure_summaries.append(summary)
    ledger = merge_ledgers(*ledgers)
    accepted, exclusions, exclusion_counts = strict_independent_filter(
        candidates, ledger, minimum_heldout_year=minimum_heldout_year
    )
    if accepted.empty:
        raise ValueError("No strictly independent expansion events passed the frozen protocol")
    accepted["expansion_event_id"] = accepted.apply(_event_identifier, axis=1)
    accepted["validation_event_id"] = accepted["expansion_event_id"]
    assert_independence(accepted, ledger, minimum_heldout_year=minimum_heldout_year)
    if accepted.expansion_event_id.duplicated().any():
        raise AssertionError("Expansion event identifiers are not unique")

    events_path = output_dir / "expansion_events.csv"
    measurements_path = output_dir / "used_measurements.csv"
    full_records_path = output_dir / "used_measurements.full.jsonl"
    provenance_path = output_dir / "natural_product_provenance.csv"
    _write_csv(accepted, events_path)
    _write_csv(exclusions, output_dir / "excluded_candidate_events.csv")
    measurements = validation_measurements(accepted, resolved).rename(
        columns={"validation_event_id": "expansion_event_id"}
    )
    _write_csv(measurements, measurements_path)
    full_count = _write_jsonl(complete_measurement_records(accepted, resolved), full_records_path)
    provenance = natural_product_provenance(accepted, resolved)
    _write_csv(provenance, provenance_path)

    source_counts = {}
    for source in ("NPASS", "COCONUT"):
        molecules = set(provenance.loc[provenance.source == source, "molecule_id"].astype(str)) if len(provenance) else set()
        source_counts[source] = {
            "events": int(accepted.molecule_id.astype(str).isin(molecules).sum()),
            "positive_events": int(accepted.loc[accepted.molecule_id.astype(str).isin(molecules), "disagreement_label"].sum()),
        }
    summary = {
        "version": EXPANSION_VERSION,
        "created_at": datetime.now(UTC).isoformat(),
        "protocol_frozen_before_outcomes": True,
        "labels_available_to_model_development": False,
        "threshold_pchembl": DISAGREEMENT_THRESHOLD,
        "candidate_events_before_independence_filter": len(candidates),
        "accepted_events": len(accepted),
        "positive_events": int(accepted.disagreement_label.sum()),
        "negative_events": int((accepted.disagreement_label == 0).sum()),
        "heldout_studies": int(accepted.heldout_study_key.nunique()),
        "targets_with_events": int(accepted.target_chembl_id.nunique()),
        "verified_natural_product_events": int(accepted.is_verified_natural_product.astype(bool).sum()),
        "verified_natural_product_positive_events": int(accepted.loc[accepted.is_verified_natural_product.astype(bool), "disagreement_label"].sum()),
        "natural_product_sources": source_counts,
        "used_measurements": len(measurements),
        "full_provenance_records": full_count,
        "exclusion_counts": exclusion_counts,
        "construction": construction,
        "previous_exposure_ledgers": exposure_summaries,
        "previous_exposure_union": {key: len(value) for key, value in ledger.items()},
        "inputs": {
            "frozen_target_specifications": {"path": str(specifications), "sha256": _sha256(specifications)},
            "resolved_target_specifications": {"path": str(resolved), "sha256": _sha256(resolved)},
        },
    }
    (output_dir / "preparation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def classify_protocol_explanation(row: pd.Series | dict[str, Any]) -> str:
    value = row if isinstance(row, dict) else row.to_dict()
    comparison_columns = sorted(
        key for key in value if key.startswith("protocol_") and key.endswith("_comparison")
    )
    statuses = [value[key] for key in comparison_columns]
    if "different" in statuses:
        return "documented_protocol_difference_possible_explanation"
    known = sum(status == "same" for status in statuses)
    both_complete = (
        value.get("prior_protocol_availability") == "complete"
        and value.get("heldout_protocol_availability") == "complete"
    )
    if both_complete and known >= 4:
        return "unexplained_with_partial_documented_similarity"
    return "unresolved_due_to_incomplete_protocol_metadata"


def _case_studies(enriched: pd.DataFrame, measurements: pd.DataFrame) -> list[dict[str, Any]]:
    cases = []
    for _, row in enriched[enriched.disagreement_label.astype(int) == 1].iterrows():
        event_measurements = measurements[measurements.expansion_event_id == row.expansion_event_id]
        def side(role: str) -> list[dict[str, Any]]:
            subset = event_measurements[event_measurements.evidence_role == role]
            return [
                {
                    "activity_id": item.activity_id,
                    "assay_chembl_id": item.assay_chembl_id,
                    "document_chembl_id": item.document_chembl_id,
                    "document_year": None if pd.isna(item.document_year) else int(item.document_year),
                    "pchembl": float(item.computed_pactivity),
                    "standard_value": float(item.standard_value),
                    "standard_units": item.standard_units,
                    "assay_description": None if pd.isna(item.assay_description) else item.assay_description,
                }
                for item in subset.itertuples(index=False)
            ]
        changed = sorted(
            column.removeprefix("protocol_").removesuffix("_comparison")
            for column in enriched.columns
            if column.startswith("protocol_") and column.endswith("_comparison") and row[column] == "different"
        )
        cases.append({
            "expansion_event_id": row.expansion_event_id,
            "target_chembl_id": row.target_chembl_id,
            "target_name": row.target_name,
            "molecule_id": row.molecule_id,
            "canonical_smiles": row.canonical_smiles,
            "verified_natural_product": bool(row.is_verified_natural_product),
            "natural_product_sources": str(row.natural_product_sources),
            "prior_median_pchembl": float(row.prior_median_pactivity),
            "independent_median_pchembl": float(row.heldout_median_pactivity),
            "signed_shift_pchembl": float(row.signed_pactivity_shift),
            "absolute_shift_pchembl": float(row.absolute_pactivity_shift),
            "threshold_pchembl": DISAGREEMENT_THRESHOLD,
            "protocol_interpretation": row.protocol_interpretation,
            "documented_differences": changed,
            "historical_measurements": side("prior"),
            "independent_measurements": side("heldout"),
            "interpretation_limit": "Association is descriptive; protocol differences do not establish causation.",
        })
    return cases


def finalize_expansion(
    event_dataset: Path,
    measurement_dataset: Path,
    raw_assays: Path,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(event_dataset)
    measurements = pd.read_csv(measurement_dataset)
    assays = load_and_extract(raw_assays)
    enriched = enrich_events(events, assays)
    enriched["protocol_interpretation"] = enriched.apply(classify_protocol_explanation, axis=1)
    _write_csv(enriched, output_dir / "enriched_expansion_events.csv")
    _write_csv(assays, output_dir / "assay_protocol_features.csv")
    cases = _case_studies(enriched, measurements)
    (output_dir / "disagreement_case_studies.json").write_text(
        json.dumps(cases, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    case_table = pd.DataFrame([{
        key: value if not isinstance(value, list) else ";".join(str(item) for item in value)
        for key, value in case.items() if key not in {"historical_measurements", "independent_measurements"}
    } for case in cases])
    _write_csv(case_table, output_dir / "disagreement_case_studies.csv")

    interpretation = Counter(enriched.loc[enriched.disagreement_label.astype(int) == 1, "protocol_interpretation"])
    comparison_coverage = {}
    for column in sorted(c for c in enriched if c.startswith("protocol_") and c.endswith("_comparison")):
        comparison_coverage[column] = dict(Counter(enriched[column].fillna("unknown")))
    coverage = {
        "events": len(enriched),
        "positive_events": int(enriched.disagreement_label.sum()),
        "negative_events": int((enriched.disagreement_label == 0).sum()),
        "heldout_studies": int(enriched.heldout_study_key.nunique()),
        "targets": int(enriched.target_chembl_id.nunique()),
        "assays_requested": len(set().union(*[
            set(str(value).split(";")) for value in pd.concat([events.prior_assay_ids, events.heldout_assay_ids])
        ])),
        "assays_retrieved": int(assays.assay_chembl_id.nunique()),
        "verified_natural_products": int(enriched.is_verified_natural_product.astype(bool).sum()),
        "verified_natural_product_positives": int(enriched.loc[enriched.is_verified_natural_product.astype(bool), "disagreement_label"].sum()),
        "positive_protocol_interpretations": dict(sorted(interpretation.items())),
        "comparison_coverage": comparison_coverage,
    }
    (output_dir / "coverage_statistics.json").write_text(json.dumps(coverage, indent=2, sort_keys=True), encoding="utf-8")
    report = [
        "# Assay Disagreement Dataset Expansion", "",
        f"Generated {datetime.now(UTC).date().isoformat()} under `{EXPANSION_VERSION}`.", "",
        "## Cohort", "",
        f"The strictly independent expansion contains **{coverage['events']} compound-target events** across **{coverage['heldout_studies']} held-out publications** and **{coverage['targets']} targets**. "
        f"There are **{coverage['positive_events']} disagreements** (absolute prior-versus-independent median shift >= 1.0 pChEMBL) and **{coverage['negative_events']} negatives**. "
        f"Exact standardized-structure matches identify **{coverage['verified_natural_products']} verified natural products**, including **{coverage['verified_natural_product_positives']} disagreements**.", "",
        "## Scientific interpretation", "",
        "All labels use exact IC50 measurements within one target and a frozen compatible assay context. Historical and held-out evidence have disjoint ChEMBL assay and document identifiers. Incompatible endpoints, censored relations, normalization mismatches, potential duplicates, and records carrying ChEMBL validity comments are excluded rather than treated as disagreements.", "",
        "Protocol comparisons are post-outcome explanatory annotations. A documented difference is a possible explanation only; it is not causal evidence. `unexplained_with_partial_documented_similarity` means at least four reported dimensions agree while other protocol fields may remain unknown. `unresolved_due_to_incomplete_protocol_metadata` means ChEMBL does not report enough protocol detail to interpret the shift. Missing fields are never inferred.", "",
        "## Positive-case protocol classification", "",
    ]
    for label, count in sorted(interpretation.items()):
        report.append(f"- `{label}`: {count}")
    report += [
        "", "## Label quarantine", "",
        "These new labels are reserved as an independent dataset. They were not used to train, tune, select, or evaluate a model in this phase. No model artifact is loaded by this pipeline.", "",
        "## Limitations", "",
        "- ChEMBL assay descriptions are often summaries, not full experimental protocols.",
        "- Structural membership in NPASS/COCONUT verifies database membership, not biosynthetic origin or clinical efficacy.",
        "- Latest-year evidence can include several records from one publication and is therefore one study-level event, not several independent confirmations.",
        "- A one-log-unit threshold is a prespecified operational definition of disagreement, not proof of biological failure.",
    ]
    (output_dir / "SCIENTIFIC_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {
        "version": EXPANSION_VERSION,
        "created_at": datetime.now(UTC).isoformat(),
        "extraction_version": EXTRACTION_VERSION,
        "no_model_training_or_evaluation": True,
        "label_quarantine": "not available to model development",
        "inputs": {path.name: _sha256(path) for path in (event_dataset, measurement_dataset, raw_assays)},
        "outputs": {},
    }
    for name in (
        "enriched_expansion_events.csv", "assay_protocol_features.csv",
        "disagreement_case_studies.json", "disagreement_case_studies.csv",
        "coverage_statistics.json", "SCIENTIFIC_REPORT.md",
    ):
        manifest["outputs"][name] = _sha256(output_dir / name)
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return coverage


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent assay-disagreement dataset expansion")
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--specs", type=Path, required=True)
    prepare.add_argument("--data-root", type=Path, required=True)
    prepare.add_argument("--previous-specs", type=Path, action="append", required=True)
    prepare.add_argument("--previous-events", type=Path, action="append", required=True)
    prepare.add_argument("--output-dir", type=Path, required=True)
    prepare.add_argument("--minimum-heldout-year", type=int, default=2020)
    finalize = commands.add_parser("finalize")
    finalize.add_argument("--events", type=Path, required=True)
    finalize.add_argument("--measurements", type=Path, required=True)
    finalize.add_argument("--assays", type=Path, required=True)
    finalize.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_expansion(
            args.specs, args.data_root, args.previous_specs, args.previous_events,
            args.output_dir, minimum_heldout_year=args.minimum_heldout_year,
        )
    else:
        result = finalize_expansion(args.events, args.measurements, args.assays, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
