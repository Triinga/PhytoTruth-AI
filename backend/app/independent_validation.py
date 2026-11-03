from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import joblib
import numpy as np
import pandas as pd

from .assay_disagreement import (
    _jsonl,
    _metric_values,
    _sha256,
    _write_csv,
    build_dataset,
    clustered_metric_intervals,
    clustered_model_comparison,
)
from .multi_target_benchmark import load_target_specs


VALIDATION_VERSION = "independent-assay-disagreement-validation-v1"


def _ids(value: Any) -> set[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    return {item.strip() for item in str(value).split(";") if item.strip()}


def build_previous_exposure_ledger(
    previous_specs: Path, previous_events: Path
) -> tuple[dict[str, set[str]], dict[str, Any]]:
    """Index every identifier in prior source snapshots, not only modeled rows."""
    ledger = {
        "targets": set(), "molecules": set(), "activities": set(),
        "assays": set(), "documents": set(),
    }
    source_counts: dict[str, Any] = {}
    for spec in load_target_specs(previous_specs):
        target = str(spec["target_chembl_id"])
        ledger["targets"].add(target)
        count = 0
        for record in _jsonl(Path(spec["activities"])):
            count += 1
            if record.get("molecule_chembl_id"):
                ledger["molecules"].add(str(record["molecule_chembl_id"]))
            activity = record.get("activity_id", record.get("source_record_id"))
            if activity is not None:
                ledger["activities"].add(str(activity))
            if record.get("assay_chembl_id"):
                ledger["assays"].add(str(record["assay_chembl_id"]))
            if record.get("document_chembl_id"):
                ledger["documents"].add(str(record["document_chembl_id"]))
        source_counts[target] = count

    # Defense in depth: include every provenance identifier persisted in the
    # prior 531-event table even if a future source format changes.
    events = pd.read_csv(previous_events)
    for value in events.get("molecule_id", pd.Series(dtype=str)):
        ledger["molecules"].add(str(value))
    for column, key in (
        ("prior_activity_ids", "activities"), ("heldout_activity_ids", "activities"),
        ("prior_assay_ids", "assays"), ("heldout_assay_ids", "assays"),
        ("prior_document_ids", "documents"), ("heldout_document_ids", "documents"),
    ):
        if column in events:
            for value in events[column]:
                ledger[key].update(_ids(value))
    summary = {
        "prior_targets": sorted(ledger["targets"]),
        "prior_source_records_by_target": source_counts,
        "unique_molecules": len(ledger["molecules"]),
        "unique_activities": len(ledger["activities"]),
        "unique_assays": len(ledger["assays"]),
        "unique_documents": len(ledger["documents"]),
        "previous_specs_sha256": _sha256(previous_specs),
        "previous_events_sha256": _sha256(previous_events),
    }
    return ledger, summary


def event_overlap_reasons(row: pd.Series | dict[str, Any], ledger: dict[str, set[str]]) -> list[str]:
    value = row if isinstance(row, dict) else row.to_dict()
    reasons = []
    if str(value.get("target_chembl_id")) in ledger["targets"]:
        reasons.append("previous_target")
    if str(value.get("molecule_id")) in ledger["molecules"]:
        reasons.append("previous_compound")
    if (_ids(value.get("prior_activity_ids")) | _ids(value.get("heldout_activity_ids"))) & ledger["activities"]:
        reasons.append("previous_activity")
    if (_ids(value.get("prior_assay_ids")) | _ids(value.get("heldout_assay_ids"))) & ledger["assays"]:
        reasons.append("previous_assay")
    if (_ids(value.get("prior_document_ids")) | _ids(value.get("heldout_document_ids"))) & ledger["documents"]:
        reasons.append("previous_publication")
    return reasons


def strict_independent_filter(
    candidates: pd.DataFrame,
    ledger: dict[str, set[str]],
    *,
    minimum_heldout_year: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    exclusions: list[dict[str, Any]] = []
    accepted_indices: list[int] = []
    duplicated_new_molecules = set(
        candidates.loc[candidates.duplicated("molecule_id", keep=False), "molecule_id"].astype(str)
    )
    for index, row in candidates.iterrows():
        reasons = event_overlap_reasons(row, ledger)
        if int(row["heldout_year"]) < minimum_heldout_year:
            reasons.append("heldout_before_frozen_temporal_boundary")
        if str(row["partition"]) != "evaluation":
            reasons.append("not_evaluation_partition")
        if str(row["molecule_id"]) in duplicated_new_molecules:
            reasons.append("compound_repeated_across_new_targets")
        if _ids(row["prior_assay_ids"]) & _ids(row["heldout_assay_ids"]):
            reasons.append("within_event_assay_overlap")
        if _ids(row["prior_document_ids"]) & _ids(row["heldout_document_ids"]):
            reasons.append("within_event_publication_overlap")
        if reasons:
            exclusions.append({
                "target_chembl_id": row["target_chembl_id"],
                "molecule_id": row["molecule_id"],
                "heldout_year": int(row["heldout_year"]),
                "reasons": ";".join(sorted(set(reasons))),
            })
        else:
            accepted_indices.append(index)
    accepted = candidates.loc[accepted_indices].copy().reset_index(drop=True)
    rejected = pd.DataFrame(exclusions)
    counts: Counter[str] = Counter()
    for value in rejected.get("reasons", pd.Series(dtype=str)):
        counts.update(str(value).split(";"))
    for reason in (
        "previous_target", "previous_compound", "previous_activity", "previous_assay",
        "previous_publication", "heldout_before_frozen_temporal_boundary",
        "not_evaluation_partition", "compound_repeated_across_new_targets",
        "within_event_assay_overlap", "within_event_publication_overlap",
    ):
        counts.setdefault(reason, 0)
    return accepted, rejected, dict(sorted(counts.items()))


def assert_independence(frame: pd.DataFrame, ledger: dict[str, set[str]], *, minimum_heldout_year: int) -> None:
    if frame.molecule_id.duplicated().any():
        raise AssertionError("Independent validation contains duplicate compounds")
    if (frame.heldout_year.astype(int) < minimum_heldout_year).any():
        raise AssertionError("Independent validation violates temporal boundary")
    for _, row in frame.iterrows():
        reasons = event_overlap_reasons(row, ledger)
        if reasons:
            raise AssertionError(f"Prior-exposure leakage detected: {reasons}")
        if _ids(row.prior_assay_ids) & _ids(row.heldout_assay_ids):
            raise AssertionError("Within-event assay leakage")
        if _ids(row.prior_document_ids) & _ids(row.heldout_document_ids):
            raise AssertionError("Within-event publication leakage")


def validation_measurements(frame: pd.DataFrame, specs: Path) -> pd.DataFrame:
    """Materialize only measurements used in accepted events with event/role provenance."""
    activity_roles: dict[str, tuple[str, str]] = {}
    for _, event in frame.iterrows():
        event_id = str(event["validation_event_id"])
        for activity in _ids(event["prior_activity_ids"]):
            activity_roles[activity] = (event_id, "prior")
        for activity in _ids(event["heldout_activity_ids"]):
            activity_roles[activity] = (event_id, "heldout")
    rows = []
    for spec in load_target_specs(specs):
        for record in _jsonl(Path(spec["activities"])):
            activity = str(record.get("activity_id", record.get("source_record_id")))
            if activity not in activity_roles:
                continue
            event_id, role = activity_roles[activity]
            raw = record.get("raw_source") or {}
            rows.append({
                "validation_event_id": event_id,
                "evidence_role": role,
                "target_chembl_id": spec["target_chembl_id"],
                "molecule_chembl_id": record.get("molecule_chembl_id"),
                "activity_id": activity,
                "assay_chembl_id": record.get("assay_chembl_id"),
                "document_chembl_id": record.get("document_chembl_id"),
                "document_year": raw.get("document_year"),
                "document_journal": raw.get("document_journal"),
                "assay_description": raw.get("assay_description"),
                "assay_type": record.get("assay_type"),
                "target_organism": record.get("target_organism"),
                "bao_format": record.get("bao_format"),
                "standard_type": record.get("standard_type"),
                "standard_relation": record.get("standard_relation"),
                "standard_value": record.get("standard_value"),
                "standard_units": record.get("standard_units"),
                "computed_pactivity": record.get("computed_pactivity"),
                "canonical_smiles": record.get("canonical_smiles"),
                "assay_context_key": record.get("assay_context_key"),
                "source": record.get("source", "ChEMBL"),
            })
    found = {row["activity_id"] for row in rows}
    missing = sorted(set(activity_roles) - found)
    if missing:
        raise AssertionError(f"Accepted event measurements missing from normalized sources: {missing[:5]}")
    return pd.DataFrame(rows).sort_values(
        ["validation_event_id", "evidence_role", "document_year", "activity_id"]
    ).reset_index(drop=True)


def _load_frozen_models(full_path: Path, selected_path: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    checksums = {"full_model": _sha256(full_path), "selected_model": _sha256(selected_path)}
    full = joblib.load(full_path)
    selected = joblib.load(selected_path)
    required_full = {"model", "features", "validation_status"}
    if not isinstance(full, dict) or not required_full <= set(full):
        raise ValueError("Frozen full-model artifact has an unsupported schema")
    if not isinstance(selected, dict) or selected.get("selected_model") != "constant_prevalence":
        raise ValueError("Frozen study-generalized artifact is not the selected constant-prevalence model")
    return full, selected, checksums


def score_frozen_models(
    frame: pd.DataFrame,
    full_artifact: dict[str, Any],
    selected_artifact: dict[str, Any],
) -> pd.DataFrame:
    predictions = frame.copy()
    features = list(full_artifact["features"])
    required = set(features) | {"target_chembl_id"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Independent dataset is missing frozen model features: {sorted(missing)}")
    predictions["probability_frozen_full"] = full_artifact["model"].predict_proba(
        frame[features + ["target_chembl_id"]]
    )[:, 1]
    predictions["probability_frozen_selected"] = float(selected_artifact["constant_prevalence"])
    predictions["frozen_full_validation_status"] = str(full_artifact["validation_status"])
    predictions["frozen_selected_validation_status"] = str(selected_artifact["validation_status"])
    return predictions


def _evaluate_scope(
    frame: pd.DataFrame, probability_column: str, *, iterations: int, seed: int
) -> dict[str, Any]:
    point = _metric_values(frame.disagreement_label, frame[probability_column])
    intervals = clustered_metric_intervals(
        frame, probability_column, iterations=iterations, seed=seed
    ) if len(frame) and frame.heldout_study_key.nunique() >= 2 else {
        name: None for name in ("auroc", "pr_auc", "brier", "log_loss", "ece_10")
    }
    return {"point": point, "study_clustered_ci95": intervals}


def evaluate_frozen_predictions(
    predictions: pd.DataFrame, *, iterations: int, seed: int
) -> dict[str, Any]:
    result: dict[str, Any] = {"overall": {}, "natural_products": {}, "targets": {}}
    natural = predictions[predictions.is_verified_natural_product.astype(bool)]
    for index, (name, column) in enumerate((
        ("frozen_full", "probability_frozen_full"),
        ("frozen_selected", "probability_frozen_selected"),
    )):
        result["overall"][name] = _evaluate_scope(
            predictions, column, iterations=iterations, seed=seed + index
        )
        result["natural_products"][name] = _evaluate_scope(
            natural, column, iterations=iterations, seed=seed + 100 + index
        ) if len(natural) else {"point": _metric_values([], []), "study_clustered_ci95": {}}
    for target, group in predictions.groupby("target_chembl_id", sort=True):
        result["targets"][str(target)] = {
            name: _evaluate_scope(group, column, iterations=iterations, seed=seed + 200 + index)
            for index, (name, column) in enumerate((
                ("frozen_full", "probability_frozen_full"),
                ("frozen_selected", "probability_frozen_selected"),
            ))
        }
    result["paired_full_minus_selected"] = clustered_model_comparison(
        predictions, "probability_frozen_full", "probability_frozen_selected",
        iterations=iterations, seed=seed + 500,
    )
    return result


def _dataset_summary(frame: pd.DataFrame) -> dict[str, Any]:
    natural = frame[frame.is_verified_natural_product.astype(bool)]
    by_target = []
    for target, group in frame.groupby("target_chembl_id", sort=True):
        group_natural = group[group.is_verified_natural_product.astype(bool)]
        by_target.append({
            "target_chembl_id": str(target),
            "target_name": str(group.target_name.iloc[0]),
            "events": len(group),
            "disagreements": int(group.disagreement_label.sum()),
            "studies": int(group.heldout_study_key.nunique()),
            "natural_products": len(group_natural),
            "natural_product_disagreements": int(group_natural.disagreement_label.sum()),
            "heldout_year_min": int(group.heldout_year.min()),
            "heldout_year_max": int(group.heldout_year.max()),
        })
    return {
        "events": len(frame),
        "disagreements": int(frame.disagreement_label.sum()),
        "prevalence": float(frame.disagreement_label.mean()) if len(frame) else None,
        "studies": int(frame.heldout_study_key.nunique()),
        "targets": int(frame.target_chembl_id.nunique()),
        "natural_products": len(natural),
        "natural_product_disagreements": int(natural.disagreement_label.sum()),
        "natural_product_studies": int(natural.heldout_study_key.nunique()),
        "by_target": by_target,
    }


def _support_decision(summary: dict[str, Any], metrics: dict[str, Any], protocol: dict[str, Any]) -> dict[str, Any]:
    minimum = protocol["minimum_overall_evaluation"]
    enough = (
        summary["events"] >= int(minimum["events"])
        and summary["disagreements"] >= int(minimum["disagreements"])
        and summary["studies"] >= int(minimum["heldout_studies"])
    )
    full = metrics["overall"]["frozen_full"]
    selected = metrics["overall"]["frozen_selected"]
    auroc_ci = full["study_clustered_ci95"].get("auroc")
    pr_ci = full["study_clustered_ci95"].get("pr_auc")
    criteria = {
        "minimum_dataset_met": enough,
        "auroc_ci_lower_above_half": bool(auroc_ci and auroc_ci[0] > 0.5),
        "pr_auc_ci_lower_above_prevalence": bool(pr_ci and pr_ci[0] > summary["prevalence"]),
        "brier_no_worse_than_frozen_selected": bool(
            full["point"]["brier"] <= selected["point"]["brier"]
        ),
    }
    supported = all(criteria.values())
    return {
        "status": "independently_supported" if supported else "not_independently_supported",
        "prespecified_rule": protocol["support_rule"],
        "criteria": criteria,
        "overall_minimum_deficits": {
            "events": max(0, int(minimum["events"]) - summary["events"]),
            "disagreements": max(0, int(minimum["disagreements"]) - summary["disagreements"]),
            "heldout_studies": max(0, int(minimum["heldout_studies"]) - summary["studies"]),
        },
        "model_or_threshold_tuning_performed": False,
    }


def _write_svg(path: Path, metrics: dict[str, Any], prevalence: float) -> None:
    full = metrics["overall"]["frozen_full"]["point"]
    selected = metrics["overall"]["frozen_selected"]["point"]
    values = [
        ("AUROC", float(full["auroc"]), float(selected["auroc"])),
        ("PR-AUC", float(full["pr_auc"]), float(selected["pr_auc"])),
        ("Brier (lower better)", float(full["brier"]), float(selected["brier"])),
    ]
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="470" viewBox="0 0 1000 470">',
        '<rect width="100%" height="100%" fill="#07131d"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#dcecf4}.muted{fill:#8da9b8}.full{fill:#2dd4bf}.base{fill:#f59e0b}</style>',
        '<text x="45" y="48" font-size="24" font-weight="700">Independent frozen-model validation</text>',
        f'<text x="45" y="75" font-size="14" class="muted">Observed disagreement prevalence: {prevalence:.3f}</text>',
    ]
    for index, (label, full_value, selected_value) in enumerate(values):
        y = 125 + index * 105
        lines.extend([
            f'<text x="45" y="{y}" font-size="16">{html.escape(label)}</text>',
            f'<rect class="full" x="250" y="{y - 22}" width="{600 * full_value:.1f}" height="25" rx="4"/>',
            f'<text x="865" y="{y - 4}" font-size="14">full {full_value:.3f}</text>',
            f'<rect class="base" x="250" y="{y + 13}" width="{600 * selected_value:.1f}" height="25" rx="4"/>',
            f'<text x="865" y="{y + 32}" font-size="14">frozen fallback {selected_value:.3f}</text>',
        ])
    lines.append('</svg>')
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_report(path: Path, results: dict[str, Any]) -> None:
    summary = results["dataset"]
    metrics = results["metrics"]
    full = metrics["overall"]["frozen_full"]
    selected = metrics["overall"]["frozen_selected"]
    natural = metrics["natural_products"]["frozen_full"]
    rows = []
    for target in summary["by_target"]:
        rows.append(
            f"| {target['target_chembl_id']} | {target['target_name']} | {target['events']} | "
            f"{target['disagreements']} | {target['studies']} | {target['natural_products']} | "
            f"{target['natural_product_disagreements']} |"
        )
    ci = full["study_clustered_ci95"]
    lines = [
        "# Independent assay-disagreement validation",
        "",
        "## Independence contract",
        "",
        f"This analyst-independent target replication used **{summary['targets']} previously unused ChEMBL targets** and held-out publications from 2020 onward. Every compound, activity, assay, and publication appearing anywhere in the previous five normalized target snapshots was excluded. The frozen model artifacts were loaded without fitting, calibration, threshold selection, or parameter changes; their SHA-256 values were verified unchanged after scoring.",
        "",
        "All validation targets are outside the frozen model's five-target training domain. The existing one-hot encoder therefore treats target identity as unknown; results measure cross-target as well as cross-study generalization.",
        "",
        "## Dataset",
        "",
        f"The final dataset contains **{summary['events']} compounds**, **{summary['disagreements']} disagreements**, and **{summary['studies']} held-out publications**. It includes **{summary['natural_products']} structure-verified natural products** with **{summary['natural_product_disagreements']} disagreements across {summary['natural_product_studies']} studies**. These events are backed by **{results['source_measurements']['records']} traceable measurements** ({results['source_measurements']['prior_records']} prior and {results['source_measurements']['heldout_records']} held out).",
        "",
        "| Target | Name | Events | Disagreements | Studies | Natural products | NP disagreements |",
        "|---|---|---:|---:|---:|---:|---:|",
        *rows,
        "",
        "## Frozen-model performance",
        "",
        "| Model | AUROC | Study-clustered 95% CI | PR-AUC | Study-clustered 95% CI | Brier | ECE-10 | Calibration intercept | Calibration slope |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        f"| Frozen full logistic | {full['point']['auroc']} | {ci.get('auroc')} | {full['point']['pr_auc']} | {ci.get('pr_auc')} | {full['point']['brier']} | {full['point']['ece_10']} | {full['point']['calibration_intercept']} | {full['point']['calibration_slope']} |",
        f"| Frozen development-selected fallback | {selected['point']['auroc']} | {selected['study_clustered_ci95'].get('auroc')} | {selected['point']['pr_auc']} | {selected['study_clustered_ci95'].get('pr_auc')} | {selected['point']['brier']} | {selected['point']['ece_10']} | {selected['point']['calibration_intercept']} | {selected['point']['calibration_slope']} |",
        "",
        "## Natural-product subgroup",
        "",
        f"The frozen full model had AUROC `{natural['point']['auroc']}`, PR-AUC `{natural['point']['pr_auc']}`, Brier `{natural['point']['brier']}`, and study-clustered intervals `{natural['study_clustered_ci95']}`. The prespecified subgroup sufficiency criteria are reported in `results.json`; metrics remain descriptive if those criteria are not met.",
        "",
        "## Decision",
        "",
        f"**{results['validation_decision']['status']}**. Criteria: `{json.dumps(results['validation_decision']['criteria'], sort_keys=True)}`.",
        "",
        f"The prespecified overall minimum was not met by `{results['validation_decision']['overall_minimum_deficits']}`. Natural-product deficits were `{results['validation_decision']['natural_product_minimum_deficits']}`. These deficits describe the minimum additional eligible observations needed for the planned analysis; they do not guarantee adequate statistical power.",
        "",
        "This result concerns prediction of later-assay disagreement only. It does not identify which assay is correct and does not imply drug efficacy, safety, biological failure, or clinical performance.",
        "",
        "## Limitations",
        "",
        "- The targets were frozen without inspecting disagreement outcomes, but chosen partly for ChEMBL activity volume and biological relevance; they are not a random sample of pharmacological targets.",
        "- All targets are unseen by the frozen target encoder, so cross-target shift is deliberately severe and target-specific intercepts cannot transfer.",
        "- Publication identifiers are imperfect proxies for independent laboratories and protocols.",
        "- Natural-product membership is structure-verified NPASS/COCONUT membership, not proof of natural origin for the specific assayed sample.",
        "- The latest-publication construction is retrospective and ChEMBL deposition can be incomplete or delayed.",
        "- Confidence intervals preserve held-out-study clustering within target; few clusters or positives can still make them unstable.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m app.pipeline independent-assay-validation --specs ../experiments/independent_validation_targets.json --protocol ../experiments/independent_validation_protocol.json --previous-specs ../experiments/multitarget_targets.json --previous-events ../experiments/assay_disagreement_generalization/dataset.csv --full-model ../experiments/assay_disagreement/assay_disagreement_model.joblib --selected-model ../experiments/assay_disagreement_generalization/study_generalized_model.joblib --output-dir ../experiments/independent_assay_validation --bootstrap-iterations 2000 --seed 2026",
        "```",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_independent_validation(
    specs: Path,
    protocol_path: Path,
    previous_specs: Path,
    previous_events: Path,
    full_model_path: Path,
    selected_model_path: Path,
    output_dir: Path,
    *,
    bootstrap_iterations: int = 2000,
    seed: int = 2026,
) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol.get("protocol_version") != VALIDATION_VERSION or not protocol.get("frozen_before_outcome_construction"):
        raise ValueError("Independent validation protocol is absent or not frozen")
    if bootstrap_iterations < 1:
        raise ValueError("bootstrap_iterations must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)

    full_artifact, selected_artifact, model_hashes_before = _load_frozen_models(
        full_model_path, selected_model_path
    )
    candidates, candidate_manifest = build_dataset(specs)
    ledger, ledger_summary = build_previous_exposure_ledger(previous_specs, previous_events)
    independent, exclusions, exclusion_counts = strict_independent_filter(
        candidates, ledger, minimum_heldout_year=int(protocol["minimum_heldout_year"])
    )
    if independent.empty:
        raise ValueError("No independent events remain after prespecified leakage exclusions")
    independent.insert(0, "validation_event_id", [
        hashlib.sha256(
            f"{row.target_chembl_id}|{row.molecule_id}|{row.heldout_study_key}".encode()
        ).hexdigest()[:20]
        for row in independent.itertuples(index=False)
    ])
    assert_independence(
        independent, ledger, minimum_heldout_year=int(protocol["minimum_heldout_year"])
    )
    measurements = validation_measurements(independent, specs)
    predictions = score_frozen_models(independent, full_artifact, selected_artifact)
    metrics = evaluate_frozen_predictions(
        predictions, iterations=bootstrap_iterations, seed=seed
    )
    summary = _dataset_summary(predictions)
    minimum_np = protocol["minimum_natural_product_subgroup"]
    np_sufficient = (
        summary["natural_products"] >= int(minimum_np["events"])
        and summary["natural_product_disagreements"] >= int(minimum_np["disagreements"])
        and summary["natural_product_studies"] >= int(minimum_np["heldout_studies"])
    )
    decision = _support_decision(summary, metrics, protocol)
    decision["natural_product_subgroup_sufficient"] = bool(np_sufficient)
    decision["natural_product_minimum_deficits"] = {
        "events": max(0, int(minimum_np["events"]) - summary["natural_products"]),
        "disagreements": max(0, int(minimum_np["disagreements"]) - summary["natural_product_disagreements"]),
        "heldout_studies": max(0, int(minimum_np["heldout_studies"]) - summary["natural_product_studies"]),
    }
    model_hashes_after = {
        "full_model": _sha256(full_model_path), "selected_model": _sha256(selected_model_path)
    }
    if model_hashes_before != model_hashes_after:
        raise AssertionError("A frozen model artifact changed during validation")

    results = {
        "generated_at": datetime.now(UTC).isoformat(),
        "validation_version": VALIDATION_VERSION,
        "seed": seed,
        "bootstrap_iterations": bootstrap_iterations,
        "dataset": summary,
        "source_measurements": {
            "records": len(measurements),
            "prior_records": int((measurements.evidence_role == "prior").sum()),
            "heldout_records": int((measurements.evidence_role == "heldout").sum()),
        },
        "candidate_events_before_independence_filter": len(candidates),
        "exclusions": exclusion_counts,
        "previous_exposure_ledger": ledger_summary,
        "metrics": metrics,
        "validation_decision": decision,
        "model_hashes_before": model_hashes_before,
        "model_hashes_after": model_hashes_after,
        "models_retrained_or_tuned": False,
        "candidate_manifest": candidate_manifest,
    }
    _write_csv(output_dir / "independent_dataset.csv", independent)
    _write_csv(output_dir / "validation_measurements.csv", measurements)
    _write_csv(output_dir / "frozen_model_predictions.csv", predictions)
    _write_csv(output_dir / "independence_exclusions.csv", exclusions)
    target_metrics_rows = []
    for target, model_results in metrics["targets"].items():
        for model, values in model_results.items():
            target_metrics_rows.append({
                "target_chembl_id": target, "model": model,
                **values["point"],
                **{f"{name}_ci95": json.dumps(interval) for name, interval in values["study_clustered_ci95"].items()},
            })
    _write_csv(output_dir / "target_metrics.csv", pd.DataFrame(target_metrics_rows))
    (output_dir / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "created_at": datetime.now(UTC).isoformat(),
        "validation_version": VALIDATION_VERSION,
        "inputs": {
            "specs": {"path": str(specs.resolve()), "sha256": _sha256(specs)},
            "protocol": {"path": str(protocol_path.resolve()), "sha256": _sha256(protocol_path)},
            "previous_specs": {"path": str(previous_specs.resolve()), "sha256": _sha256(previous_specs)},
            "previous_events": {"path": str(previous_events.resolve()), "sha256": _sha256(previous_events)},
            "full_model": {"path": str(full_model_path.resolve()), "sha256": model_hashes_before["full_model"]},
            "selected_model": {"path": str(selected_model_path.resolve()), "sha256": model_hashes_before["selected_model"]},
        },
        "new_target_inputs": [
            {
                "target_chembl_id": spec["target_chembl_id"],
                "activities": spec["activities"], "activities_sha256": _sha256(Path(spec["activities"])),
                "natural_products": spec["natural_products"], "natural_products_sha256": _sha256(Path(spec["natural_products"])),
            }
            for spec in load_target_specs(specs)
        ],
        "outputs": {},
    }
    _write_svg(output_dir / "independent_validation.svg", metrics, float(summary["prevalence"]))
    _write_report(output_dir / "SCIENTIFIC_REPORT.md", results)
    for path in (
        output_dir / "independent_dataset.csv", output_dir / "validation_measurements.csv",
        output_dir / "frozen_model_predictions.csv",
        output_dir / "independence_exclusions.csv", output_dir / "target_metrics.csv",
        output_dir / "results.json", output_dir / "independent_validation.svg",
        output_dir / "SCIENTIFIC_REPORT.md",
    ):
        manifest["outputs"][path.name] = _sha256(path)
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    command = (
        f'python -m app.pipeline independent-assay-validation --specs "{specs.resolve()}" '
        f'--protocol "{protocol_path.resolve()}" --previous-specs "{previous_specs.resolve()}" '
        f'--previous-events "{previous_events.resolve()}" --full-model "{full_model_path.resolve()}" '
        f'--selected-model "{selected_model_path.resolve()}" --output-dir "{output_dir.resolve()}" '
        f'--bootstrap-iterations {bootstrap_iterations} --seed {seed}'
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent frozen-model assay-disagreement validation")
    parser.add_argument("--specs", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--previous-specs", type=Path, required=True)
    parser.add_argument("--previous-events", type=Path, required=True)
    parser.add_argument("--full-model", type=Path, required=True)
    parser.add_argument("--selected-model", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    results = run_independent_validation(
        args.specs, args.protocol, args.previous_specs, args.previous_events,
        args.full_model, args.selected_model, args.output_dir,
        bootstrap_iterations=args.bootstrap_iterations, seed=args.seed,
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "dataset": results["dataset"],
        "validation_decision": results["validation_decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
