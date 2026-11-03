from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from scipy.stats import norm
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from .reversal import METRIC_NAMES, evaluate_ranking
from .reversal_experiment import (
    _jsonl,
    _sha256,
    build_benchmark_rows,
    latest_publication_split,
    natural_product_membership,
    select_compatible_measurements,
)


PRIMARY_METRICS = ("ndcg_at_k", "kendall_tau")
BASELINE_METHODS = (
    "median_potency",
    "baseline_max_potency",
    "mean_potency",
    "latest_training_potency",
    "lower_confidence_bound",
    "phytotruth_full",
)
FEATURE_SETS: dict[str, tuple[str, ...]] = {
    "minimal": (
        "mean_minus_median",
        "latest_minus_median",
        "log_documents",
        "mad",
    ),
    "reliability": (
        "log_measurements",
        "log_assays",
        "log_documents",
        "mad",
        "pactivity_range",
        "contradiction_rate",
        "duplicate_rate",
        "reliability_heuristic",
        "confidence_width",
    ),
    "full": (
        "mean_minus_median",
        "latest_minus_median",
        "log_measurements",
        "log_assays",
        "log_documents",
        "mad",
        "pactivity_range",
        "contradiction_rate",
        "duplicate_rate",
        "reliability_heuristic",
        "confidence_width",
        "history_year_span",
    ),
}
ABLATION_REMOVALS: dict[str, tuple[str, ...]] = {
    "without_support": ("log_measurements", "log_assays", "log_documents"),
    "without_consistency": ("mad", "pactivity_range", "contradiction_rate"),
    "without_reliability": ("duplicate_rate", "reliability_heuristic", "confidence_width"),
    "without_recent_drift": ("mean_minus_median", "latest_minus_median", "history_year_span"),
}


def discover_targets(data_roots: Sequence[Path]) -> list[dict[str, Any]]:
    """Discover every local target with both normalized and NP crossmatch files."""
    specs: list[dict[str, Any]] = []
    seen: set[str] = set()
    for root in sorted(Path(path).resolve() for path in data_roots):
        for activities in sorted(root.glob("CHEMBL*.normalized.jsonl")):
            target_id = activities.name.split(".", 1)[0]
            natural = activities.with_name(f"{target_id}.natural.jsonl")
            if not natural.is_file() or target_id in seen:
                continue
            first = next(_jsonl(activities), None)
            if not first or str(first.get("target_chembl_id")) != target_id:
                raise ValueError(f"Target identity mismatch in {activities}")
            raw = first.get("raw_source") or {}
            specs.append({
                "target_chembl_id": target_id,
                "target_name": raw.get("target_pref_name") or target_id,
                "activities": activities,
                "natural_products": natural,
            })
            seen.add(target_id)
    if not specs:
        raise ValueError("No paired target datasets were discovered")
    return sorted(specs, key=lambda item: item["target_chembl_id"])


def deterministic_target_partition(
    target_ids: Sequence[str], *, development_fraction: float, seed: int
) -> tuple[list[str], list[str]]:
    """Partition whole biological targets without looking at outcomes or scores."""
    unique = sorted(set(target_ids))
    if len(unique) != len(target_ids):
        raise ValueError("target IDs must be unique")
    if not 0 < development_fraction < 1 or len(unique) < 3:
        raise ValueError("a valid fraction and at least three targets are required")
    ordered = sorted(
        unique,
        key=lambda target: hashlib.sha256(f"{seed}:{target}".encode()).hexdigest(),
    )
    n_development = min(len(ordered) - 1, max(1, round(len(ordered) * development_fraction)))
    return sorted(ordered[:n_development]), sorted(ordered[n_development:])


def _features(row: dict[str, Any]) -> dict[str, float]:
    return {
        "mean_minus_median": float(row["mean_potency"]) - float(row["median_potency"]),
        "latest_minus_median": float(row["latest_training_potency"]) - float(row["median_potency"]),
        "log_measurements": math.log1p(float(row["train_measurements"])),
        "log_assays": math.log1p(float(row["train_assays"])),
        "log_documents": math.log1p(float(row["train_documents"])),
        "mad": float(row["mad"]),
        "pactivity_range": float(row["pactivity_range"]),
        "contradiction_rate": float(row["contradiction_rate"]),
        "duplicate_rate": float(row["duplicate_rate"]),
        "reliability_heuristic": float(row["reliability_heuristic"]),
        "confidence_width": float(row["confidence_interval_high"]) - float(row["confidence_interval_low"]),
        "history_year_span": float(row["train_last_year"]) - float(row["train_first_year"]),
    }


def _matrix(rows: Sequence[dict[str, Any]], feature_names: Sequence[str]) -> np.ndarray:
    return np.asarray([[ _features(row)[name] for name in feature_names] for row in rows], dtype=float)


def _fit_ridge(
    target_rows: dict[str, Sequence[dict[str, Any]]], feature_names: Sequence[str], alpha: float
) -> tuple[StandardScaler, Ridge]:
    rows = [row for target in sorted(target_rows) for row in target_rows[target]]
    weights = [1.0 / len(target_rows[target]) for target in sorted(target_rows) for _ in target_rows[target]]
    x = _matrix(rows, feature_names)
    y = np.asarray([
        float(row["heldout_median_pactivity"]) - float(row["median_potency"])
        for row in rows
    ])
    scaler = StandardScaler().fit(x)
    model = Ridge(alpha=alpha).fit(scaler.transform(x), y, sample_weight=np.asarray(weights))
    return scaler, model


def _predict(
    rows: Sequence[dict[str, Any]], feature_names: Sequence[str], scaler: StandardScaler, model: Ridge
) -> list[float]:
    residual = model.predict(scaler.transform(_matrix(rows, feature_names)))
    return [float(row["median_potency"]) + float(delta) for row, delta in zip(rows, residual)]


def _target_metrics(
    target_rows: dict[str, Sequence[dict[str, Any]]], method: str, *, k: int, active_threshold: float
) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for target, rows in sorted(target_rows.items()):
        output[target] = asdict(evaluate_ranking(
            [float(row[method]) for row in rows],
            [float(row["heldout_median_pactivity"]) for row in rows],
            min(k, len(rows)),
            active_threshold=active_threshold,
        ))
    return output


def _exact_sign_flip_p(differences: Sequence[float]) -> float:
    values = np.asarray(differences, dtype=float)
    if not len(values):
        return 1.0
    observed = abs(float(np.mean(values)))
    if observed <= 1e-15:
        return 1.0
    if len(values) <= 20:
        exceedances = 0
        total = 2 ** len(values)
        for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
            exceedances += int(abs(float(np.mean(values * np.asarray(signs)))) >= observed - 1e-15)
        return exceedances / total
    rng = np.random.default_rng(2026)
    iterations = 200_000
    signs = rng.choice((-1.0, 1.0), size=(iterations, len(values)))
    exceedances = int(np.count_nonzero(np.abs(np.mean(signs * values, axis=1)) >= observed - 1e-15))
    return (exceedances + 1) / (iterations + 1)


def _holm_adjust(p_values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(p_values, key=p_values.get)
    adjusted: dict[str, float] = {}
    running = 0.0
    count = len(ordered)
    for index, name in enumerate(ordered):
        running = max(running, (count - index) * p_values[name])
        adjusted[name] = min(1.0, running)
    return adjusted


def target_cluster_inference(
    target_rows: dict[str, Sequence[dict[str, Any]]],
    candidate: str,
    *,
    baseline_method: str = "median_potency",
    k: int,
    active_threshold: float,
    bootstrap_iterations: int,
    seed: int,
) -> dict[str, Any]:
    """Inference at the target level; molecules are not treated as independent targets."""
    baseline = _target_metrics(target_rows, baseline_method, k=k, active_threshold=active_threshold)
    contender = _target_metrics(target_rows, candidate, k=k, active_threshold=active_threshold)
    targets = sorted(target_rows)
    rng = np.random.default_rng(seed)
    raw_p: dict[str, float] = {}
    result: dict[str, Any] = {}
    for metric in METRIC_NAMES:
        left = np.asarray([baseline[target][metric] for target in targets])
        right = np.asarray([contender[target][metric] for target in targets])
        differences = right - left
        bootstrap = [
            float(np.mean(differences[rng.integers(0, len(targets), size=len(targets))]))
            for _ in range(bootstrap_iterations)
        ]
        point = float(np.mean(differences))
        p_value = _exact_sign_flip_p(differences)
        raw_p[metric] = p_value
        sd = float(np.std(differences, ddof=1)) if len(differences) > 1 else 0.0
        dz = point / sd if sd > 0 else (math.inf if point > 0 else -math.inf if point < 0 else 0.0)
        hedges = dz * (1 - 3 / (4 * len(differences) - 5)) if len(differences) > 1 and math.isfinite(dz) else dz
        result[metric] = {
            "baseline": float(np.mean(left)),
            "candidate": float(np.mean(right)),
            "difference": point,
            "difference_ci95": [float(value) for value in np.quantile(bootstrap, [0.025, 0.975])],
            "target_sign_flip_p": p_value,
            "target_wins": int(np.count_nonzero(differences > 0)),
            "target_ties": int(np.count_nonzero(np.isclose(differences, 0))),
            "target_losses": int(np.count_nonzero(differences < 0)),
            "hedges_g_paired": hedges,
            "target_differences": {target: float(delta) for target, delta in zip(targets, differences)},
        }
    adjusted = _holm_adjust(raw_p)
    for metric in METRIC_NAMES:
        result[metric]["holm_adjusted_p"] = adjusted[metric]
    return result


def _attach_predictions(
    target_rows: dict[str, Sequence[dict[str, Any]]],
    method: str,
    feature_names: Sequence[str],
    scaler: StandardScaler,
    model: Ridge,
) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = {}
    for target, rows in target_rows.items():
        values = _predict(rows, feature_names, scaler, model)
        output[target] = [{**row, method: value} for row, value in zip(rows, values)]
    return output


def cross_validate_candidates(
    development: dict[str, Sequence[dict[str, Any]]],
    *,
    alphas: Sequence[float],
    k: int,
    active_threshold: float,
    bootstrap_iterations: int,
    seed: int,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    grid: list[dict[str, Any]] = []
    predictions_by_config: dict[str, dict[str, list[dict[str, Any]]]] = {}
    target_ids = sorted(development)
    for feature_set, feature_names in FEATURE_SETS.items():
        for alpha in alphas:
            config = f"{feature_set}:ridge:{alpha:g}"
            predictions: dict[str, list[dict[str, Any]]] = {}
            for heldout_target in target_ids:
                training = {target: rows for target, rows in development.items() if target != heldout_target}
                scaler, model = _fit_ridge(training, feature_names, alpha)
                rows = development[heldout_target]
                values = _predict(rows, feature_names, scaler, model)
                predictions[heldout_target] = [
                    {**row, "optimized_ranking": value} for row, value in zip(rows, values)
                ]
            inference = target_cluster_inference(
                predictions,
                "optimized_ranking",
                k=k,
                active_threshold=active_threshold,
                bootstrap_iterations=bootstrap_iterations,
                seed=seed + len(grid),
            )
            objective = (
                inference["ndcg_at_k"]["difference"]
                + 0.25 * inference["kendall_tau"]["difference"]
            )
            grid.append({
                "config": config,
                "feature_set": feature_set,
                "features": list(feature_names),
                "alpha": alpha,
                "development_objective": objective,
                "development_inference": inference,
            })
            predictions_by_config[config] = predictions
    selected = max(grid, key=lambda row: (row["development_objective"], -row["alpha"], row["config"]))
    inference = selected["development_inference"]
    passed = all(
        inference[metric]["difference_ci95"][0] > 0
        and inference[metric]["holm_adjusted_p"] < 0.05
        for metric in PRIMARY_METRICS
    )
    selection = {
        "candidate_config": selected["config"],
        "feature_set": selected["feature_set"],
        "features": selected["features"],
        "alpha": selected["alpha"],
        "passed_development_gate": passed,
        "deployed_method": "optimized_ranking" if passed else "median_potency",
        "gate": "Both target-macro nDCG@10 and Kendall tau must have target-bootstrap CI lower bounds above zero and Holm-adjusted exact target sign-flip p<0.05 in leave-one-development-target-out validation.",
        "development_inference": inference,
        "grid": grid,
    }
    return selection, predictions_by_config[selected["config"]]


def _power_summary(inference: dict[str, Any], n_targets: int, alpha: float = 0.05) -> dict[str, Any]:
    primary = inference["ndcg_at_k"]
    effect = float(primary["hedges_g_paired"])
    required = None
    if math.isfinite(effect) and effect > 0:
        required = math.ceil(((norm.ppf(1 - alpha / 2) + norm.ppf(0.8)) / effect) ** 2)
    return {
        "independent_evaluation_targets": n_targets,
        "minimum_attainable_two_sided_exact_p": 2 / (2 ** n_targets),
        "observed_paired_standardized_effect": effect,
        "approximate_targets_for_80_percent_power_at_alpha_0_05": required,
        "note": "The sample-size approximation uses a paired-normal approximation and is planning guidance, not a post-hoc proof of power.",
    }


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    materialized = list(rows)
    if not materialized:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(materialized[0]))
        writer.writeheader()
        writer.writerows(materialized)


def _report(results: dict[str, Any]) -> str:
    primary = results["evaluation_inference"]["ndcg_at_k"]
    tau = results["evaluation_inference"]["kendall_tau"]
    status = results["scientific_conclusion"]["status"]
    selection = results["development_selection"]
    power = results["power_analysis"]
    return f"""# Robust multi-target PhytoTruth evaluation

Generated {results['generated_at']} from real ChEMBL IC50 measurements already present in the repository.

## Leakage contract

Whole biological targets were deterministically assigned to development or evaluation before any outcome-based model fitting. For every molecule, the latest publication year was held out and all scoring features came from strictly earlier, document-disjoint evidence. Ridge strength and feature family were selected only by leave-one-development-target-out validation. The evaluation targets were used once for the frozen comparison.

This is a retrospective database benchmark, not a prospective or independent-laboratory replication. Several source targets have appeared in earlier repository analyses, so the result is not described as a new confirmatory study.

## Coverage

- Discovered targets: {results['coverage']['discovered_targets']}
- Eligible targets: {results['coverage']['eligible_targets']}
- Development targets: {results['coverage']['development_targets']}
- Evaluation targets: {results['coverage']['evaluation_targets']}
- Evaluation molecules with future-publication evidence: {results['coverage']['evaluation_molecules']}

## Development-only selection

The best nested-development candidate was `{selection['candidate_config']}`. It {'passed' if selection['passed_development_gate'] else 'did not pass'} the frozen two-metric gate, so the deployed ranking is `{selection['deployed_method']}`.

## Locked evaluation

Optimized candidate minus median potency macro nDCG@10: {primary['difference']:.4f} (95% target-bootstrap CI {primary['difference_ci95'][0]:.4f} to {primary['difference_ci95'][1]:.4f}; exact target sign-flip p={primary['target_sign_flip_p']:.6f}; Holm-adjusted p={primary['holm_adjusted_p']:.6f}; Hedges paired g={primary['hedges_g_paired']:.3f}; wins/ties/losses={primary['target_wins']}/{primary['target_ties']}/{primary['target_losses']}).

Kendall tau difference: {tau['difference']:.4f} (95% CI {tau['difference_ci95'][0]:.4f} to {tau['difference_ci95'][1]:.4f}; exact p={tau['target_sign_flip_p']:.6f}; Holm-adjusted p={tau['holm_adjusted_p']:.6f}).

Scientific status: **{status}**. {results['scientific_conclusion']['summary']}

## Statistical power

With {power['independent_evaluation_targets']} target clusters, the minimum attainable two-sided exact sign-flip p-value is {power['minimum_attainable_two_sided_exact_p']:.6f}. The observed paired standardized effect was {power['observed_paired_standardized_effect']:.3f}. Approximate target count for 80% power at alpha 0.05: {power['approximate_targets_for_80_percent_power_at_alpha_0_05']}.

## Ablations

`ablation_inference.json` and `ablation_summary.csv` compare the frozen full candidate with removal of support, consistency, reliability/uncertainty, and recent-drift features. These are secondary analyses; Holm-adjusted target-level p-values are reported and negative contributions are retained.

## Interpretation limits

- ChEMBL publications are proxies for independent experiments and may share laboratories or copied measurements.
- Target-level inference avoids treating molecules as independent biological targets, but available targets are still a convenience sample.
- nDCG@10 measures future-publication ranking agreement, not efficacy, safety, or clinical utility.
- Natural-product membership is never a training feature; subgroup counts are descriptive because target-held-out NP power remains limited.
"""


def run_evaluation(
    data_roots: Sequence[Path],
    output_dir: Path,
    *,
    development_fraction: float = 0.6,
    minimum_molecules: int = 25,
    k: int = 10,
    active_threshold: float = 6.0,
    bootstrap_iterations: int = 10_000,
    seed: int = 2026,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    specs = discover_targets(data_roots)
    prepared: dict[str, list[dict[str, Any]]] = {}
    target_summary: list[dict[str, Any]] = []
    input_manifest: list[dict[str, Any]] = []
    for spec in specs:
        records = list(_jsonl(spec["activities"]))
        selected, context, exclusions = select_compatible_measurements(records)
        splits, split_exclusions = latest_publication_split(selected)
        rows = build_benchmark_rows(splits, natural_product_membership(spec["natural_products"]))
        eligible = len(rows) >= minimum_molecules
        target_id = spec["target_chembl_id"]
        if eligible:
            prepared[target_id] = rows
        target_summary.append({
            "target_chembl_id": target_id,
            "target_name": spec["target_name"],
            "eligible": eligible,
            "normalized_measurements": len(records),
            "compatible_measurements": len(selected),
            "future_evidence_molecules": len(rows),
            "verified_natural_products": sum(bool(row["is_verified_natural_product"]) for row in rows),
            "assay_context_key": context,
            "coverage_exclusions": json.dumps(exclusions, sort_keys=True),
            "split_exclusions": json.dumps(split_exclusions, sort_keys=True),
        })
        input_manifest.append({
            "target_chembl_id": target_id,
            "activities": str(spec["activities"]),
            "activities_sha256": _sha256(spec["activities"]),
            "natural_products": str(spec["natural_products"]),
            "natural_products_sha256": _sha256(spec["natural_products"]),
        })
    if len(prepared) < 6:
        raise ValueError(f"Only {len(prepared)} targets met coverage; at least six are required")
    development_ids, evaluation_ids = deterministic_target_partition(
        list(prepared), development_fraction=development_fraction, seed=seed
    )
    development = {target: prepared[target] for target in development_ids}
    evaluation = {target: prepared[target] for target in evaluation_ids}
    alphas = (0.01, 0.1, 1.0, 10.0, 100.0)
    selection, development_predictions = cross_validate_candidates(
        development,
        alphas=alphas,
        k=k,
        active_threshold=active_threshold,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed,
    )
    feature_names = tuple(selection["features"])
    scaler, model = _fit_ridge(development, feature_names, float(selection["alpha"]))
    evaluation_predictions = _attach_predictions(
        evaluation, "optimized_ranking", feature_names, scaler, model
    )
    evaluation_inference = target_cluster_inference(
        evaluation_predictions,
        "optimized_ranking",
        k=k,
        active_threshold=active_threshold,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed + 50_000,
    )
    current_inference = target_cluster_inference(
        evaluation,
        "phytotruth_full",
        k=k,
        active_threshold=active_threshold,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed + 60_000,
    )

    ablations: dict[str, Any] = {}
    ablation_rows: list[dict[str, Any]] = []
    for index, (name, removed) in enumerate(ABLATION_REMOVALS.items()):
        removed_from_selected = tuple(feature for feature in removed if feature in feature_names)
        remaining = tuple(feature for feature in feature_names if feature not in removed_from_selected)
        if not removed_from_selected:
            ablations[name] = {
                "removed_features": [],
                "status": "not_applicable_to_selected_feature_set",
                "full_minus_ablation": None,
            }
            continue
        if not remaining:
            continue
        ablation_scaler, ablation_model = _fit_ridge(development, remaining, float(selection["alpha"]))
        ablation_predictions = _attach_predictions(evaluation, name, remaining, ablation_scaler, ablation_model)
        combined = {
            target: [
                {**full, name: ablated[name]}
                for full, ablated in zip(evaluation_predictions[target], ablation_predictions[target])
            ]
            for target in evaluation_ids
        }
        inference = target_cluster_inference(
            combined,
            "optimized_ranking",
            baseline_method=name,
            k=k,
            active_threshold=active_threshold,
            bootstrap_iterations=bootstrap_iterations,
            seed=seed + 70_000 + index,
        )
        ablations[name] = {
            "removed_features": list(removed_from_selected),
            "status": "executed",
            "full_minus_ablation": inference,
        }
        for metric in METRIC_NAMES:
            ablation_rows.append({
                "ablation": name,
                "metric": metric,
                "full_minus_ablation": inference[metric]["difference"],
                "ci95_low": inference[metric]["difference_ci95"][0],
                "ci95_high": inference[metric]["difference_ci95"][1],
                "exact_p": inference[metric]["target_sign_flip_p"],
                "holm_adjusted_p": inference[metric]["holm_adjusted_p"],
            })

    primary = evaluation_inference["ndcg_at_k"]
    tau = evaluation_inference["kendall_tau"]
    superiority = (
        selection["passed_development_gate"]
        and primary["difference_ci95"][0] > 0
        and tau["difference_ci95"][0] > 0
        and primary["holm_adjusted_p"] < 0.05
        and tau["holm_adjusted_p"] < 0.05
    )
    status = "target_heldout_superiority_supported" if superiority else "superiority_not_demonstrated"
    summary = (
        "The optimized evidence-aware ranker passed development selection and independently improved both prespecified ranking metrics."
        if superiority
        else "The frozen target-held-out evaluation does not support a superiority claim; the median-potency guardrail remains deployed."
    )
    results: dict[str, Any] = {
        "experiment": "Robust target-held-out PhytoTruth evaluation",
        "generated_at": datetime.now(UTC).isoformat(),
        "seed": seed,
        "k": k,
        "active_threshold": active_threshold,
        "bootstrap_iterations": bootstrap_iterations,
        "independence_unit": "biological_target",
        "coverage": {
            "discovered_targets": len(specs),
            "eligible_targets": len(prepared),
            "development_targets": len(development),
            "evaluation_targets": len(evaluation),
            "development_molecules": sum(map(len, development.values())),
            "evaluation_molecules": sum(map(len, evaluation.values())),
            "development_target_ids": development_ids,
            "evaluation_target_ids": evaluation_ids,
        },
        "development_selection": selection,
        "evaluation_inference": evaluation_inference,
        "current_heuristic_evaluation_inference": current_inference,
        "power_analysis": _power_summary(evaluation_inference, len(evaluation)),
        "ablations": ablations,
        "scientific_conclusion": {
            "status": status,
            "superiority_claim": superiority,
            "deployed_ranking": selection["deployed_method"] if superiority else "median_potency",
            "summary": summary,
            "confirmatory_status": "retrospective_reused_database_targets_not_prospective_confirmation",
        },
    }
    target_metric_rows: list[dict[str, Any]] = []
    for metric in METRIC_NAMES:
        metric_result = evaluation_inference[metric]
        for target, difference in metric_result["target_differences"].items():
            target_metric_rows.append({"target_chembl_id": target, "metric": metric, "optimized_minus_median": difference})
    coefficient_rows = [
        {"feature": feature, "standardized_coefficient": float(value)}
        for feature, value in zip(feature_names, model.coef_)
    ]
    partitions = {target: "development" for target in development_ids} | {
        target: "evaluation" for target in evaluation_ids
    }
    for row in target_summary:
        row["partition"] = partitions.get(row["target_chembl_id"], "ineligible")
    _write_csv(output_dir / "target_coverage.csv", target_summary)
    _write_csv(output_dir / "target_metric_differences.csv", target_metric_rows)
    _write_csv(output_dir / "ablation_summary.csv", ablation_rows)
    _write_csv(output_dir / "model_coefficients.csv", coefficient_rows)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    (output_dir / "ablation_inference.json").write_text(json.dumps(ablations, indent=2, sort_keys=True), encoding="utf-8")
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    command = (
        "python -m app.robust_multitarget_evaluation "
        + " ".join(f"--data-root {path}" for path in data_roots)
        + f" --output-dir {output_dir} --development-fraction {development_fraction}"
        + f" --minimum-molecules {minimum_molecules} --k {k}"
        + f" --active-threshold {active_threshold:g} --bootstrap-iterations {bootstrap_iterations} --seed {seed}\n"
    )
    (output_dir / "command.txt").write_text(command, encoding="utf-8")
    manifest = {
        "inputs": input_manifest,
        "results_sha256": _sha256(output_dir / "results.json"),
        "target_partition_uses_outcomes": False,
        "model_selection_data": "development targets only",
        "evaluation_access": "single frozen evaluation after nested development selection",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Run robust target-held-out PhytoTruth evaluation")
    parser.add_argument("--data-root", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-fraction", type=float, default=0.6)
    parser.add_argument("--minimum-molecules", type=int, default=25)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--active-threshold", type=float, default=6.0)
    parser.add_argument("--bootstrap-iterations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    results = run_evaluation(
        args.data_root,
        args.output_dir,
        development_fraction=args.development_fraction,
        minimum_molecules=args.minimum_molecules,
        k=args.k,
        active_threshold=args.active_threshold,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    print(json.dumps({
        "status": results["scientific_conclusion"]["status"],
        "coverage": results["coverage"],
        "ndcg": results["evaluation_inference"]["ndcg_at_k"],
    }, indent=2))


if __name__ == "__main__":
    main()
