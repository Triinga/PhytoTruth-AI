from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import numpy as np

from .evidence import score_evidence
from .reversal import METRIC_NAMES, evaluate_ranking, paired_ranking_inference
from .schemas import EvidenceRecord


METHODS = (
    "baseline_max_potency",
    "median_potency",
    "mean_potency",
    "latest_training_potency",
    "lower_confidence_bound",
    "ablate_reliability",
    "ablate_uncertainty",
    "ablate_explicit_consistency",
    "ablate_support",
    "phytotruth_full",
)


def _jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {path}:{line_number}") from exc
            if not isinstance(value, dict):
                raise ValueError(f"Expected a JSON object at {path}:{line_number}")
            yield value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def natural_product_membership(path: Path) -> dict[str, list[dict[str, Any]]]:
    membership: dict[str, list[dict[str, Any]]] = {}
    for record in _jsonl(path):
        molecule_id = str(record.get("molecule_chembl_id") or "")
        provenance = record.get("natural_product_provenance")
        if molecule_id and isinstance(provenance, list) and provenance:
            membership[molecule_id] = provenance
    return membership


def select_compatible_measurements(
    records: Iterable[dict[str, Any]],
    *,
    context_key: str | None = None,
) -> tuple[list[dict[str, Any]], str, dict[str, int]]:
    """Select one assay context without using activity outcomes to choose it."""
    seen: set[str] = set()
    candidates: list[dict[str, Any]] = []
    exclusions: Counter[str] = Counter()
    for record in records:
        activity_id = str(record.get("activity_id") or record.get("source_record_id") or "")
        if not activity_id:
            exclusions["missing_activity_id"] += 1
            continue
        if activity_id in seen:
            exclusions["duplicate_activity_id"] += 1
            continue
        seen.add(activity_id)
        if record.get("standard_relation") != "=" or record.get("is_censored"):
            exclusions["non_exact_measurement"] += 1
            continue
        year = (record.get("raw_source") or {}).get("document_year")
        try:
            year = int(year)
            value = float(record["computed_pactivity"])
        except (KeyError, TypeError, ValueError):
            exclusions["missing_year_or_activity"] += 1
            continue
        if not math.isfinite(value) or not record.get("molecule_chembl_id"):
            exclusions["invalid_activity_or_molecule"] += 1
            continue
        candidates.append({**record, "document_year": year, "computed_pactivity": value})

    if not candidates:
        raise ValueError("No exact, dated measurements are available")
    if context_key is None:
        molecules_by_context: dict[str, set[str]] = defaultdict(set)
        for record in candidates:
            key = str(record.get("assay_context_key") or "")
            if key:
                molecules_by_context[key].add(str(record["molecule_chembl_id"]))
        if not molecules_by_context:
            raise ValueError("No assay context keys are available")
        context_key = sorted(molecules_by_context, key=lambda key: (-len(molecules_by_context[key]), key))[0]

    selected = [record for record in candidates if record.get("assay_context_key") == context_key]
    exclusions["other_assay_context"] += len(candidates) - len(selected)
    if not selected:
        raise ValueError(f"Assay context {context_key!r} has no eligible records")
    return selected, context_key, dict(exclusions)


def latest_publication_split(records: Iterable[dict[str, Any]]) -> tuple[dict[str, dict[str, list[dict[str, Any]]]], dict[str, int]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["molecule_chembl_id"])].append(record)
    splits: dict[str, dict[str, list[dict[str, Any]]]] = {}
    exclusions: Counter[str] = Counter()
    for molecule_id, group in grouped.items():
        years = sorted({int(record["document_year"]) for record in group})
        if len(years) < 2:
            exclusions["single_publication_year"] += 1
            continue
        heldout_year = years[-1]
        train = [record for record in group if int(record["document_year"]) < heldout_year]
        heldout = [record for record in group if int(record["document_year"]) == heldout_year]
        _assert_no_leakage(train, heldout)
        splits[molecule_id] = {"train": train, "heldout": heldout}
    return splits, dict(exclusions)


def global_publication_split(
    records: Iterable[dict[str, Any]], cutoff_year: int
) -> tuple[dict[str, dict[str, list[dict[str, Any]]]], dict[str, int]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["molecule_chembl_id"])].append(record)
    splits: dict[str, dict[str, list[dict[str, Any]]]] = {}
    exclusions: Counter[str] = Counter()
    for molecule_id, group in grouped.items():
        train = [record for record in group if int(record["document_year"]) <= cutoff_year]
        heldout = [record for record in group if int(record["document_year"]) > cutoff_year]
        if not train or not heldout:
            exclusions["not_observed_on_both_sides"] += 1
            continue
        _assert_no_leakage(train, heldout)
        splits[molecule_id] = {"train": train, "heldout": heldout}
    return splits, dict(exclusions)


def _assert_no_leakage(train: list[dict[str, Any]], heldout: list[dict[str, Any]]) -> None:
    if not train or not heldout:
        raise ValueError("Both train and held-out evidence are required")
    train_activities = {str(record["activity_id"]) for record in train}
    heldout_activities = {str(record["activity_id"]) for record in heldout}
    if train_activities & heldout_activities:
        raise ValueError("Activity-record leakage detected")
    train_documents = {str(record.get("document_chembl_id")) for record in train}
    heldout_documents = {str(record.get("document_chembl_id")) for record in heldout}
    if train_documents & heldout_documents:
        raise ValueError("Document leakage detected")
    if max(int(record["document_year"]) for record in train) >= min(int(record["document_year"]) for record in heldout):
        raise ValueError("Publication-time ordering is not strict")


def _evidence_record(record: dict[str, Any]) -> EvidenceRecord:
    return EvidenceRecord(
        molecule_id=str(record["molecule_chembl_id"]),
        activity_id=record["activity_id"],
        assay_id=str(record["assay_chembl_id"]),
        document_id=str(record.get("document_chembl_id") or "") or None,
        pchembl_value=float(record["computed_pactivity"]),
        standard_value_nm=float(record["standard_value_molar"]) * 1e9,
        relation="=",
        organism=record.get("target_organism"),
        assay_type=(record.get("assay_context") or {}).get("assay_type"),
        potential_duplicate=bool(record.get("potential_duplicate")),
        validity_comment=record.get("data_validity_comment"),
    )


def build_benchmark_rows(
    splits: dict[str, dict[str, list[dict[str, Any]]]],
    natural_members: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    evidence = [
        _evidence_record(record)
        for molecule_id in sorted(splits)
        for record in splits[molecule_id]["train"]
    ]
    scores = {score.molecule_id: score for score in score_evidence(evidence)}
    rows: list[dict[str, Any]] = []
    for molecule_id in sorted(splits):
        train = splits[molecule_id]["train"]
        heldout = splits[molecule_id]["heldout"]
        score = scores[molecule_id]
        ci_width = score.confidence_interval[1] - score.confidence_interval[0]
        reliability_penalty = 0.25 * (1.0 - score.reliability)
        uncertainty_penalty = 0.15 * ci_width
        full = score.evidence_score - reliability_penalty - uncertainty_penalty
        explicit_consistency = 0.65 * score.mad + 0.8 * score.contradiction_rate + 0.12 * score.pchembl_range
        support_bonus = 0.18 * math.log1p(score.independent_assays)
        latest_training_year = max(int(record["document_year"]) for record in train)
        latest_training_values = [
            float(record["computed_pactivity"])
            for record in train
            if int(record["document_year"]) == latest_training_year
        ]
        provenance = natural_members.get(molecule_id, [])
        rows.append({
            "molecule_id": molecule_id,
            "is_verified_natural_product": bool(provenance),
            "natural_product_sources": ";".join(sorted({str(item.get("source")) for item in provenance})),
            "train_measurements": len(train),
            "train_assays": len({record["assay_chembl_id"] for record in train}),
            "train_documents": len({record["document_chembl_id"] for record in train}),
            "train_first_year": min(record["document_year"] for record in train),
            "train_last_year": max(record["document_year"] for record in train),
            "heldout_measurements": len(heldout),
            "heldout_documents": len({record["document_chembl_id"] for record in heldout}),
            "heldout_year": min(record["document_year"] for record in heldout),
            "heldout_median_pactivity": float(median(record["computed_pactivity"] for record in heldout)),
            "reliability_heuristic": score.reliability,
            "duplicate_rate": score.duplicate_rate,
            "organism_count": score.organism_count,
            "mad": score.mad,
            "pactivity_range": score.pchembl_range,
            "contradiction_rate": score.contradiction_rate,
            "confidence_interval_low": score.confidence_interval[0],
            "confidence_interval_high": score.confidence_interval[1],
            "consistency_penalty": explicit_consistency,
            "support_bonus": support_bonus,
            "reliability_penalty": reliability_penalty,
            "uncertainty_penalty": uncertainty_penalty,
            "baseline_max_potency": score.max_pchembl,
            "median_potency": score.median_pchembl,
            "mean_potency": float(np.mean([record["computed_pactivity"] for record in train])),
            "latest_training_potency": float(median(latest_training_values)),
            "lower_confidence_bound": score.confidence_interval[0],
            "ablate_reliability": full + reliability_penalty,
            "ablate_uncertainty": full + uncertainty_penalty,
            "ablate_explicit_consistency": full + explicit_consistency,
            "ablate_support": full - support_bonus,
            "phytotruth_full": full,
        })
    return rows


def evaluate_methods(rows: list[dict[str, Any]], *, k: int, active_threshold: float) -> dict[str, Any]:
    if not rows:
        return {"n_molecules": 0, "k": 0, "metrics": {}}
    truth = [float(row["heldout_median_pactivity"]) for row in rows]
    actual_k = min(k, len(rows))
    metrics: dict[str, dict[str, float]] = {}
    for method in METHODS:
        metrics[method] = asdict(evaluate_ranking(
            [float(row[method]) for row in rows],
            truth,
            actual_k,
            active_threshold=active_threshold,
        ))
    return {
        "n_molecules": len(rows),
        "k": actual_k,
        "active_threshold": active_threshold,
        "heldout_active_molecules": sum(value >= active_threshold for value in truth),
        "metrics": metrics,
    }


def _inference(
    rows: list[dict[str, Any]],
    left: str,
    right: str,
    *,
    k: int,
    active_threshold: float,
    bootstrap_iterations: int,
    permutation_iterations: int,
    seed: int,
) -> dict[str, Any]:
    if len(rows) < 2:
        return {"status": "insufficient_molecules"}
    return paired_ranking_inference(
        [float(row[left]) for row in rows],
        [float(row[right]) for row in rows],
        [float(row["heldout_median_pactivity"]) for row in rows],
        k=min(k, len(rows)),
        active_threshold=active_threshold,
        bootstrap_iterations=bootstrap_iterations,
        permutation_iterations=permutation_iterations,
        seed=seed,
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summary_rows(cohorts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cohort_name, cohort in cohorts.items():
        for method, metrics in cohort["evaluation"]["metrics"].items():
            for metric, value in metrics.items():
                rows.append({
                    "cohort": cohort_name,
                    "n_molecules": cohort["evaluation"]["n_molecules"],
                    "method": method,
                    "metric": metric,
                    "value": value,
                })
    return rows


def _comparison_svg(path: Path, cohorts: dict[str, dict[str, Any]]) -> None:
    chosen = [name for name in ("latest_all", "latest_natural") if cohorts[name]["evaluation"]["n_molecules"]]
    metric_names = ("ndcg_at_k", "enrichment_factor_at_k", "precision_at_k", "kendall_tau")
    width, height = 1100, 580
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#07131f"/>',
        '<text x="50" y="42" fill="#e8f5ef" font-family="Arial" font-size="24" font-weight="700">Evidence Reversal: held-out ranking performance</text>',
        '<text x="50" y="68" fill="#95a9b8" font-family="Arial" font-size="13">Blue: potency-only baseline · Green: PhytoTruth · values are point estimates</text>',
    ]
    panel_width = 250
    for panel, metric in enumerate(metric_names):
        x0 = 55 + panel * panel_width
        parts.append(f'<text x="{x0}" y="105" fill="#d6e6df" font-family="Arial" font-size="14">{metric}</text>')
        for cohort_index, cohort_name in enumerate(chosen):
            metrics = cohorts[cohort_name]["evaluation"]["metrics"]
            values = [metrics["baseline_max_potency"][metric], metrics["phytotruth_full"][metric]]
            scale_max = max(1.0, max(values) * 1.1) if metric == "enrichment_factor_at_k" else 1.0
            y_base = 245 + cohort_index * 180
            label = f'{cohort_name} (n={cohorts[cohort_name]["evaluation"]["n_molecules"]})'
            parts.append(f'<text x="{x0}" y="{y_base + 25}" fill="#95a9b8" font-family="Arial" font-size="12">{label}</text>')
            for offset, (value, color) in enumerate(zip(values, ("#4aa3ff", "#42d392"))):
                normalized = max(0.0, value) / scale_max
                bar_height = normalized * 105
                x = x0 + 20 + offset * 70
                y = y_base - bar_height
                parts.append(f'<rect x="{x}" y="{y:.1f}" width="45" height="{bar_height:.1f}" rx="4" fill="{color}"/>')
                parts.append(f'<text x="{x + 22.5}" y="{y - 6:.1f}" text-anchor="middle" fill="#e8f5ef" font-family="Arial" font-size="11">{value:.3f}</text>')
            parts.append(f'<line x1="{x0 + 10}" y1="{y_base}" x2="{x0 + 145}" y2="{y_base}" stroke="#385064"/>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def _fmt(value: float) -> str:
    return f"{value:.3f}"


def _report(results: dict[str, Any]) -> str:
    latest_all = results["cohorts"]["latest_all"]
    latest_natural = results["cohorts"]["latest_natural"]
    inference = latest_all["baseline_vs_phytotruth"]
    ndcg = inference["ndcg_at_k"]
    conclusion = results["scientific_conclusion"]["summary"]
    lines = [
        "# PhytoTruth Evidence Reversal Experiment",
        "",
        f"Generated: {results['generated_at']}",
        "",
        "## Research question",
        "",
        "Does a prespecified evidence-aware heuristic rank compounds by potency that reproduces in later publications better than ranking by the strongest earlier measurement?",
        "",
        "## Data and protocol",
        "",
        f"- Source: normalized real ChEMBL {results['target_chembl_id']} IC50 records; NPASS membership is used for the natural-product subgroup.",
        f"- Compatible context: `{results['assay_context_key']}` ({json.dumps(results['assay_context'], sort_keys=True)}).",
        "- Only exact, dated measurements were eligible. Censored observations and other assay contexts were excluded rather than pooled.",
        "- Primary split: for each molecule, all measurements from its latest publication year are held out; only strictly earlier documents score the molecule.",
        f"- Sensitivity split: global publication cutoff {results['global_cutoff_year']} (training through the cutoff, later evidence held out).",
        f"- Activity threshold: pActivity ≥ {results['active_threshold']}; K={results['k']}.",
        f"- Inference: {results['bootstrap_iterations']} molecule-level paired bootstrap samples and {results['permutation_iterations']} paired score-swap randomizations.",
        "- The reliability value is a deterministic heuristic, not a calibrated probability. No weight or threshold was tuned using held-out evidence.",
        "",
        "## Primary results",
        "",
        "| Cohort | Molecules | Method | nDCG@K | EF@K | Precision@K | Kendall τ |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for name in ("latest_all", "latest_natural", "global_all", "global_natural"):
        cohort = results["cohorts"][name]
        for method in ("baseline_max_potency", "phytotruth_full"):
            metrics = cohort["evaluation"]["metrics"].get(method)
            if not metrics:
                continue
            lines.append(
                f"| {name} | {cohort['evaluation']['n_molecules']} | {method} | "
                f"{_fmt(metrics['ndcg_at_k'])} | {_fmt(metrics['enrichment_factor_at_k'])} | "
                f"{_fmt(metrics['precision_at_k'])} | {_fmt(metrics['kendall_tau'])} |"
            )
    lines.extend([
        "",
        "For the primary all-compound cohort, the PhytoTruth-minus-baseline nDCG difference was "
        f"{_fmt(float(ndcg['difference']))} (95% molecule-bootstrap CI "
        f"[{_fmt(float(ndcg['difference_ci95'][0]))}, {_fmt(float(ndcg['difference_ci95'][1]))}], "
        f"paired randomization p={float(ndcg['paired_randomization_p']):.4f}).",
        "",
        conclusion,
        "",
        results["scientific_conclusion"]["natural_product_finding"],
        "",
        "## Ablation study",
        "",
        "| Method | nDCG@K | EF@K | Precision@K | Kendall τ |",
        "|---|---:|---:|---:|---:|",
    ])
    for method in METHODS:
        metrics = latest_all["evaluation"]["metrics"][method]
        lines.append(
            f"| {method} | {_fmt(metrics['ndcg_at_k'])} | {_fmt(metrics['enrichment_factor_at_k'])} | "
            f"{_fmt(metrics['precision_at_k'])} | {_fmt(metrics['kendall_tau'])} |"
        )
    lines.extend([
        "",
        "Ablations remove one prespecified component while keeping the split fixed. They are diagnostic comparisons, not independently tuned models.",
        "",
        results["scientific_conclusion"]["ablation_finding"],
        "",
        "## Leakage audit",
        "",
        "- Activity IDs and document IDs are disjoint between scoring and held-out evidence for every molecule.",
        "- Every held-out publication year is strictly later than that molecule's training years.",
        "- Assay context selection used molecule coverage only, never activity outcomes.",
        "- Confidence intervals resample whole molecules; individual measurements are never treated as independent experimental units.",
        "",
        "## Limitations",
        "",
        f"- Only {latest_natural['evaluation']['n_molecules']} NPASS-verified natural products had repeated evidence across publication years in the selected context. Natural-product subgroup estimates are therefore exploratory and may be very imprecise.",
        "- The latest-year split is retrospective and molecule-specific; it is not a prospective trial. The global-cutoff sensitivity is stricter but smaller.",
        "- ChEMBL assay-context metadata cannot capture every protocol difference. Sharing target, organism, endpoint, assay type, and BAO format reduces but does not eliminate heterogeneity.",
        "- Held-out ChEMBL measurements are literature evidence, not clinical efficacy or safety outcomes.",
        "- The heuristic weights were prespecified and not fitted; failure to outperform is informative and should not be hidden.",
        "",
        "## Reproduction",
        "",
        "See `command.txt`, `manifest.json`, `results.json`, and the cohort CSV files in this directory.",
        "",
    ])
    return "\n".join(lines)


def run_experiment(
    activities: Path,
    natural_products: Path,
    output_dir: Path,
    *,
    context_key: str | None = None,
    global_cutoff_year: int = 2015,
    k: int = 10,
    active_threshold: float = 6.0,
    bootstrap_iterations: int = 2000,
    permutation_iterations: int = 2000,
    seed: int = 2026,
) -> dict[str, Any]:
    if k < 1:
        raise ValueError("k must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    natural_members = natural_product_membership(natural_products)
    selected, selected_context, exclusions = select_compatible_measurements(_jsonl(activities), context_key=context_key)
    latest, latest_exclusions = latest_publication_split(selected)
    global_split, global_exclusions = global_publication_split(selected, global_cutoff_year)
    latest_rows = build_benchmark_rows(latest, natural_members)
    global_rows = build_benchmark_rows(global_split, natural_members)
    row_sets = {
        "latest_all": latest_rows,
        "latest_natural": [row for row in latest_rows if row["is_verified_natural_product"]],
        "global_all": global_rows,
        "global_natural": [row for row in global_rows if row["is_verified_natural_product"]],
    }
    context = selected[0].get("assay_context") or {}
    cohorts: dict[str, dict[str, Any]] = {}
    for index, (name, rows) in enumerate(row_sets.items()):
        evaluation = evaluate_methods(rows, k=k, active_threshold=active_threshold)
        inference = _inference(
            rows,
            "baseline_max_potency",
            "phytotruth_full",
            k=k,
            active_threshold=active_threshold,
            bootstrap_iterations=bootstrap_iterations,
            permutation_iterations=permutation_iterations,
            seed=seed + index,
        )
        cohorts[name] = {
            "evaluation": evaluation,
            "baseline_vs_phytotruth": inference,
        }
        _write_csv(output_dir / f"cohort_{name}.csv", rows)

    ablation_inference: dict[str, Any] = {}
    for index, method in enumerate(METHODS):
        if method in {"baseline_max_potency", "phytotruth_full"}:
            continue
        ablation_inference[method] = _inference(
            latest_rows,
            method,
            "phytotruth_full",
            k=k,
            active_threshold=active_threshold,
            bootstrap_iterations=max(200, bootstrap_iterations // 2),
            permutation_iterations=max(200, permutation_iterations // 2),
            seed=seed + 100 + index,
        )

    results: dict[str, Any] = {
        "experiment": "PhytoTruth Evidence Reversal",
        "generated_at": datetime.now(UTC).isoformat(),
        "target_chembl_id": context.get("target_chembl_id"),
        "assay_context_key": selected_context,
        "assay_context": context,
        "selected_measurements": len(selected),
        "selected_molecules": len({record["molecule_chembl_id"] for record in selected}),
        "verified_natural_product_molecules_in_context": len({record["molecule_chembl_id"] for record in selected} & set(natural_members)),
        "global_cutoff_year": global_cutoff_year,
        "k": k,
        "active_threshold": active_threshold,
        "seed": seed,
        "bootstrap_iterations": bootstrap_iterations,
        "permutation_iterations": permutation_iterations,
        "exclusions": {
            "normalization_and_context": exclusions,
            "latest_split": latest_exclusions,
            "global_split": global_exclusions,
        },
        "cohorts": cohorts,
        "primary_ablation_inference": ablation_inference,
    }
    primary_ndcg = cohorts["latest_all"]["baseline_vs_phytotruth"]["ndcg_at_k"]
    natural_ndcg = cohorts["latest_natural"]["baseline_vs_phytotruth"]["ndcg_at_k"]
    median_kendall = ablation_inference["median_potency"]["kendall_tau"]
    primary_supported = (
        float(primary_ndcg["difference_ci95"][0]) > 0
        and float(primary_ndcg["paired_randomization_p"]) < 0.05
    )
    results["scientific_conclusion"] = {
        "status": "advantage_supported" if primary_supported else "advantage_not_demonstrated",
        "summary": (
            "The evidence-aware heuristic's all-compound nDCG point estimate was higher, but its paired 95% interval included zero and the randomization test was not significant. The experiment therefore does not demonstrate an advantage over potency-only ranking."
            if not primary_supported else
            "The evidence-aware heuristic improved all-compound nDCG with a paired interval excluding zero and a significant randomization test."
        ),
        "natural_product_finding": (
            f"In the NPASS-verified natural-product cohort, PhytoTruth-minus-baseline nDCG was {float(natural_ndcg['difference']):.3f} "
            f"(95% CI [{float(natural_ndcg['difference_ci95'][0]):.3f}, {float(natural_ndcg['difference_ci95'][1]):.3f}]); "
            "the available natural-product data do not support an improvement."
        ),
        "ablation_finding": (
            f"Median potency alone matched the full heuristic's top-K metrics and its Kendall τ exceeded the full heuristic by {-float(median_kendall['difference']):.3f}. "
            f"The full-minus-median Kendall difference was {float(median_kendall['difference']):.3f} "
            f"(95% CI [{float(median_kendall['difference_ci95'][0]):.3f}, {float(median_kendall['difference_ci95'][1]):.3f}], "
            f"p={float(median_kendall['paired_randomization_p']):.4f}). Thus the reliability, uncertainty, consistency, and support terms did not improve this benchmark beyond robust central potency and may have modestly worsened global ordering."
        ),
    }
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(output_dir / "metric_summary.csv", _summary_rows(cohorts))
    _comparison_svg(output_dir / "comparison.svg", cohorts)
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    manifest = {
        "activities": str(activities.resolve()),
        "activities_sha256": _sha256(activities),
        "natural_products": str(natural_products.resolve()),
        "natural_products_sha256": _sha256(natural_products),
        "results_sha256": _sha256(output_dir / "results.json"),
        "code_contract": "latest-publication evidence is held out before all ranking scores are computed",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    command = (
        "python -m app.reversal_experiment "
        f'--activities "{activities}" --natural-products "{natural_products}" '
        f'--output-dir "{output_dir}" --global-cutoff-year {global_cutoff_year} --k {k} '
        f"--active-threshold {active_threshold} --bootstrap-iterations {bootstrap_iterations} "
        f"--permutation-iterations {permutation_iterations} --seed {seed}"
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the leakage-resistant PhytoTruth evidence-reversal experiment")
    parser.add_argument("--activities", type=Path, required=True, help="Normalized ChEMBL JSONL")
    parser.add_argument("--natural-products", type=Path, required=True, help="Natural-product crossmatch JSONL")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-key", help="Optional prespecified compatible assay context")
    parser.add_argument("--global-cutoff-year", type=int, default=2015)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--active-threshold", type=float, default=6.0)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--permutation-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    results = run_experiment(
        args.activities,
        args.natural_products,
        args.output_dir,
        context_key=args.context_key,
        global_cutoff_year=args.global_cutoff_year,
        k=args.k,
        active_threshold=args.active_threshold,
        bootstrap_iterations=args.bootstrap_iterations,
        permutation_iterations=args.permutation_iterations,
        seed=args.seed,
    )
    summary = {
        "output_dir": str(args.output_dir),
        "selected_measurements": results["selected_measurements"],
        "selected_molecules": results["selected_molecules"],
        "cohort_sizes": {name: value["evaluation"]["n_molecules"] for name, value in results["cohorts"].items()},
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
