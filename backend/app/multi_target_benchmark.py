from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import asdict
from datetime import UTC, datetime
from html import escape
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Sequence

import numpy as np
from scipy.stats import rankdata

from .reliability_investigation import (
    ALL_METHODS,
    component_predictive_value,
    compound_error_analysis,
    deterministic_partition,
)
from .reversal import METRIC_NAMES, evaluate_ranking
from .reversal_experiment import (
    _jsonl,
    _sha256,
    build_benchmark_rows,
    latest_publication_split,
    natural_product_membership,
    select_compatible_measurements,
)


PRIMARY_ALTERNATIVES = ("phytotruth_full", "lower_confidence_bound")
MIN_LATEST_MOLECULES = 30
MIN_LATEST_NATURAL_PRODUCTS = 5


def load_target_specs(path: Path) -> list[dict[str, Any]]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, list) or not value:
        raise ValueError("target specification must be a non-empty JSON array")
    required = {"target_chembl_id", "target_name", "activities", "natural_products"}
    seen: set[str] = set()
    specs: list[dict[str, Any]] = []
    for item in value:
        if not isinstance(item, dict) or required - set(item):
            raise ValueError(f"each target requires {sorted(required)}")
        target_id = str(item["target_chembl_id"])
        if target_id in seen:
            raise ValueError(f"duplicate target specification: {target_id}")
        seen.add(target_id)
        spec = dict(item)
        spec["activities"] = str(Path(spec["activities"]).resolve())
        spec["natural_products"] = str(Path(spec["natural_products"]).resolve())
        specs.append(spec)
    return specs


def _metric_table(rows: Sequence[dict[str, Any]], *, k: int, active_threshold: float) -> dict[str, Any]:
    truth = [float(row["heldout_median_pactivity"]) for row in rows]
    actual_k = min(k, len(rows))
    return {
        "n_molecules": len(rows),
        "natural_products": sum(bool(row["is_verified_natural_product"]) for row in rows),
        "active_molecules": sum(value >= active_threshold for value in truth),
        "k": actual_k,
        "metrics": {
            method: asdict(evaluate_ranking(
                [float(row[method]) for row in rows],
                truth,
                actual_k,
                active_threshold=active_threshold,
            ))
            for method in ALL_METHODS
        },
    }


def _macro_metrics(
    target_rows: dict[str, Sequence[dict[str, Any]]],
    method: str,
    *,
    k: int,
    active_threshold: float,
) -> dict[str, float]:
    values: dict[str, list[float]] = defaultdict(list)
    for rows in target_rows.values():
        truth = [float(row["heldout_median_pactivity"]) for row in rows]
        metrics = asdict(evaluate_ranking(
            [float(row[method]) for row in rows],
            truth,
            min(k, len(rows)),
            active_threshold=active_threshold,
        ))
        for name, value in metrics.items():
            values[name].append(float(value))
    return {name: float(np.mean(metric_values)) for name, metric_values in values.items()}


def hierarchical_macro_inference(
    target_rows: dict[str, Sequence[dict[str, Any]]],
    candidate: str,
    *,
    k: int,
    active_threshold: float,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    """Bootstrap targets then molecules; permutation swaps paired scores within target."""
    if not target_rows:
        raise ValueError("at least one target is required")
    target_ids = sorted(target_rows)
    baseline_point = _macro_metrics(target_rows, "median_potency", k=k, active_threshold=active_threshold)
    candidate_point = _macro_metrics(target_rows, candidate, k=k, active_threshold=active_threshold)
    observed = {name: candidate_point[name] - baseline_point[name] for name in METRIC_NAMES}
    bootstrap = {name: [] for name in METRIC_NAMES}
    exceedances = {name: 0 for name in METRIC_NAMES}
    rng = np.random.default_rng(seed)

    for _ in range(iterations):
        sampled_targets = rng.choice(target_ids, size=len(target_ids), replace=True)
        target_differences: dict[str, list[float]] = defaultdict(list)
        for target_id in sampled_targets:
            rows = list(target_rows[str(target_id)])
            indices = rng.integers(0, len(rows), size=len(rows))
            sampled = [rows[int(index)] for index in indices]
            truth = [float(row["heldout_median_pactivity"]) for row in sampled]
            left = asdict(evaluate_ranking(
                [float(row["median_potency"]) for row in sampled], truth, min(k, len(sampled)),
                active_threshold=active_threshold,
            ))
            right = asdict(evaluate_ranking(
                [float(row[candidate]) for row in sampled], truth, min(k, len(sampled)),
                active_threshold=active_threshold,
            ))
            for name in METRIC_NAMES:
                target_differences[name].append(float(right[name]) - float(left[name]))
        for name in METRIC_NAMES:
            bootstrap[name].append(float(np.mean(target_differences[name])))

        permuted_differences: dict[str, list[float]] = defaultdict(list)
        for target_id in target_ids:
            rows = list(target_rows[target_id])
            truth = [float(row["heldout_median_pactivity"]) for row in rows]
            left_scores = np.asarray([float(row["median_potency"]) for row in rows])
            right_scores = np.asarray([float(row[candidate]) for row in rows])
            swap = rng.random(len(rows)) < 0.5
            permuted_left = np.where(swap, right_scores, left_scores)
            permuted_right = np.where(swap, left_scores, right_scores)
            left = asdict(evaluate_ranking(permuted_left, truth, min(k, len(rows)), active_threshold=active_threshold))
            right = asdict(evaluate_ranking(permuted_right, truth, min(k, len(rows)), active_threshold=active_threshold))
            for name in METRIC_NAMES:
                permuted_differences[name].append(float(right[name]) - float(left[name]))
        for name in METRIC_NAMES:
            delta = float(np.mean(permuted_differences[name]))
            exceedances[name] += int(abs(delta) >= abs(observed[name]) - 1e-12)

    return {
        name: {
            "baseline": baseline_point[name],
            "candidate": candidate_point[name],
            "difference": observed[name],
            "difference_ci95": [float(value) for value in np.quantile(bootstrap[name], [0.025, 0.975])],
            "paired_randomization_p": (exceedances[name] + 1) / (iterations + 1),
        }
        for name in METRIC_NAMES
    }


def select_on_development(comparisons: dict[str, dict[str, Any]]) -> dict[str, Any]:
    passing: list[str] = []
    for method in PRIMARY_ALTERNATIVES:
        result = comparisons[method]
        if all(
            float(result[metric]["difference_ci95"][0]) > 0
            and float(result[metric]["paired_randomization_p"]) < 0.05
            for metric in ("ndcg_at_k", "kendall_tau")
        ):
            passing.append(method)
    selected = max(
        passing,
        key=lambda method: comparisons[method]["ndcg_at_k"]["difference"],
        default="median_potency",
    )
    return {
        "selected_method": selected,
        "passing_candidates": passing,
        "gate": (
            "Across-target macro nDCG@K and Kendall tau-b must both have paired hierarchical-bootstrap "
            "95% interval lower bounds above zero and paired score-swap p<0.05 on development targets/molecules."
        ),
    }


def _spearman(left: Sequence[float], right: Sequence[float]) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return 0.0
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def assay_transferability(
    records: Sequence[dict[str, Any]],
    molecule_ids: set[str],
    *,
    minimum_shared_molecules: int = 10,
) -> list[dict[str, Any]]:
    """Compare contexts without ever pooling their activity measurements."""
    context_values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    context_metadata: dict[str, dict[str, Any]] = {}
    for record in records:
        molecule_id = str(record["molecule_chembl_id"])
        if molecule_id not in molecule_ids:
            continue
        key = str(record["assay_context_key"])
        context_values[key][molecule_id].append(float(record["computed_pactivity"]))
        context_metadata[key] = dict(record.get("assay_context") or {})

    rows: list[dict[str, Any]] = []
    context_keys = sorted(context_values)
    for left_index, left_key in enumerate(context_keys):
        for right_key in context_keys[left_index + 1:]:
            shared = sorted(set(context_values[left_key]) & set(context_values[right_key]))
            if len(shared) < minimum_shared_molecules:
                continue
            left = [median(context_values[left_key][molecule]) for molecule in shared]
            right = [median(context_values[right_key][molecule]) for molecule in shared]
            left_meta = context_metadata[left_key]
            right_meta = context_metadata[right_key]
            rows.append({
                "left_context": left_key,
                "right_context": right_key,
                "shared_molecules": len(shared),
                "same_organism": left_meta.get("target_organism") == right_meta.get("target_organism"),
                "same_assay_type": left_meta.get("assay_type") == right_meta.get("assay_type"),
                "same_bao_format": left_meta.get("bao_format") == right_meta.get("bao_format"),
                "left_organism": left_meta.get("target_organism"),
                "right_organism": right_meta.get("target_organism"),
                "left_assay_type": left_meta.get("assay_type"),
                "right_assay_type": right_meta.get("assay_type"),
                "left_bao_format": left_meta.get("bao_format"),
                "right_bao_format": right_meta.get("bao_format"),
                "spearman_rho": _spearman(left, right),
                "median_absolute_delta": float(np.median(np.abs(np.asarray(left) - np.asarray(right)))),
                "threshold_contradiction_rate": float(np.mean(
                    (np.asarray(left) >= 6.0) != (np.asarray(right) >= 6.0)
                )),
            })
    return rows


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _flatten_metrics(results: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for target_id, target in results["targets"].items():
        if target["status"] != "included":
            continue
        for partition in ("development", "evaluation"):
            for subgroup in ("all", "natural_products"):
                table = target["partitions"][partition][subgroup]
                for method, metrics in table["metrics"].items():
                    for metric, value in metrics.items():
                        output.append({
                            "target_chembl_id": target_id,
                            "target_name": target["target_name"],
                            "partition": partition,
                            "subgroup": subgroup,
                            "n_molecules": table["n_molecules"],
                            "method": method,
                            "metric": metric,
                            "value": value,
                        })
    return output


def _comparison_svg(path: Path, results: dict[str, Any]) -> None:
    included = [(target_id, target) for target_id, target in results["targets"].items() if target["status"] == "included"]
    methods = ("median_potency", "phytotruth_full", "lower_confidence_bound")
    labels = ("Median", "PhytoTruth", "Lower CI")
    colours = ("#12b76a", "#f04438", "#6172f3")
    width = 980
    height = 120 + 105 * len(included)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="#fff"/>', '<text x="30" y="32" font-family="sans-serif" font-size="18" fill="#101828">Untouched evaluation nDCG@10 by target</text>']
    for row_index, (target_id, target) in enumerate(included):
        y = 60 + row_index * 105
        parts.append(f'<text x="30" y="{y + 17}" font-family="sans-serif" font-size="13" fill="#344054">{escape(target_id)} · {escape(target["target_name"])}</text>')
        for method_index, method in enumerate(methods):
            value = float(target["partitions"]["evaluation"]["all"]["metrics"][method]["ndcg_at_k"])
            x = 300
            bar_y = y + method_index * 24
            parts.append(f'<rect x="{x}" y="{bar_y}" width="{value * 580:.1f}" height="16" rx="3" fill="{colours[method_index]}"/>')
            parts.append(f'<text x="{x + value * 580 + 7:.1f}" y="{bar_y + 12}" font-family="sans-serif" font-size="11" fill="#475467">{labels[method_index]} {value:.3f}</text>')
    parts.append('</svg>')
    path.write_text("".join(parts), encoding="utf-8")


def _report(results: dict[str, Any]) -> str:
    included = [target for target in results["targets"].values() if target["status"] == "included"]
    target_lines = []
    for target in included:
        evaluation = target["partitions"]["evaluation"]["all"]
        natural = target["partitions"]["evaluation"]["natural_products"]
        target_lines.append(
            f"- {target['target_chembl_id']} ({target['target_name']}): {evaluation['n_molecules']} evaluation molecules "
            f"({natural['n_molecules']} natural products); median/current/lower-CI nDCG@{evaluation['k']} = "
            f"{evaluation['metrics']['median_potency']['ndcg_at_k']:.3f}/"
            f"{evaluation['metrics']['phytotruth_full']['ndcg_at_k']:.3f}/"
            f"{evaluation['metrics']['lower_confidence_bound']['ndcg_at_k']:.3f}."
        )
    macro = results["evaluation_macro_comparisons"]
    current = macro["phytotruth_full"]
    lower = macro["lower_confidence_bound"]
    return f"""# PhytoTruth multi-target evidence benchmark

Generated {results['generated_at']} from real normalized ChEMBL IC50 measurements and structure-verified NPASS/COCONUT membership.

## Design

Five human single-protein targets were prespecified before outcomes were inspected. A target was included only if its largest compatible target/organism/endpoint/assay-type/BAO context had at least {MIN_LATEST_MOLECULES} molecules with strictly earlier, document-disjoint training publications and at least {MIN_LATEST_NATURAL_PRODUCTS} structure-verified natural products. Target selection used coverage, never activity values or ranking results.

Each included target used a latest-publication temporal holdout followed by a deterministic molecule-ID development/evaluation split stratified by natural-product membership. Methods were selected only from development results. Evaluation results were opened after the development gate was frozen. The prior median-potency guardrail remained the default unless an evidence-aware method improved both macro nDCG and Kendall tau-b with intervals excluding zero and p<0.05.

## Coverage

{chr(10).join(target_lines)}

The combined natural-product index contains {results['natural_product_index']['unique_structures']} unique standardized structures with {results['natural_product_index']['npass_records']} NPASS and {results['natural_product_index']['coconut_records']} COCONUT provenance records. Cross-source duplicate structures retain both provenance records.

## Independent evaluation

The development gate selected `{results['development_selection']['selected_method']}`. Across targets, current PhytoTruth minus median nDCG was {current['ndcg_at_k']['difference']:.3f} (95% CI {current['ndcg_at_k']['difference_ci95'][0]:.3f} to {current['ndcg_at_k']['difference_ci95'][1]:.3f}; p={current['ndcg_at_k']['paired_randomization_p']:.4f}), and Kendall difference was {current['kendall_tau']['difference']:.3f}. Lower-confidence-bound minus median nDCG was {lower['ndcg_at_k']['difference']:.3f} (95% CI {lower['ndcg_at_k']['difference_ci95'][0]:.3f} to {lower['ndcg_at_k']['difference_ci95'][1]:.3f}; p={lower['ndcg_at_k']['paired_randomization_p']:.4f}). These results do not support superiority unless the recorded development gate and independent evaluation both do so.

## Transferability

`assay_transferability.csv` compares shared compounds across incompatible contexts without pooling their measurements. It reports rank correlation, median absolute pActivity shift, and activity-threshold contradictions, stratified by organism, assay type, and BAO format agreement. These comparisons are descriptive and do not assume that evidence transfers across contexts.

Organism transfer was not estimable: the prespecified targets and all retained target contexts were human. ChEMBL activity records used here expose target organism but not a consistently populated experimental-system organism field, so no cross-organism claim is made.

## Limitations

- Evaluation is literature-temporal and molecule-partitioned, but not a prospective or external laboratory validation.
- Publication IDs do not prove laboratory independence; protocols, batches, and replicate structure remain incompletely recorded.
- Natural-product subgroup power varies by target despite combining NPASS and COCONUT.
- Target-level inference has only {len(included)} biological targets, so hierarchical intervals remain imprecise.
- Held-out bioactivity is not evidence of clinical efficacy or safety.
"""


def run_benchmark(
    specs_path: Path,
    natural_product_database: Path,
    output_dir: Path,
    *,
    development_fraction: float = 0.6,
    k: int = 10,
    active_threshold: float = 6.0,
    iterations: int = 1000,
    seed: int = 2026,
) -> dict[str, Any]:
    import sqlite3

    output_dir.mkdir(parents=True, exist_ok=True)
    specs = load_target_specs(specs_path)
    targets: dict[str, Any] = {}
    partition_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
    all_transfer_rows: list[dict[str, Any]] = []
    all_error_rows: list[dict[str, Any]] = []
    all_harmful_pairs: list[dict[str, Any]] = []
    all_component_rows: list[dict[str, Any]] = []
    all_cohort_rows: list[dict[str, Any]] = []
    dataset_rows: list[dict[str, Any]] = []
    input_manifest: list[dict[str, Any]] = []

    for target_index, spec in enumerate(specs):
        target_id = str(spec["target_chembl_id"])
        activities_path = Path(spec["activities"])
        natural_path = Path(spec["natural_products"])
        natural_members = natural_product_membership(natural_path)
        all_records = list(_jsonl(activities_path))
        selected, context_key, context_exclusions = select_compatible_measurements(all_records)
        splits, split_exclusions = latest_publication_split(selected)
        rows = build_benchmark_rows(splits, natural_members)
        natural_count = sum(bool(row["is_verified_natural_product"]) for row in rows)
        latest_source_counts: dict[str, int] = defaultdict(int)
        for row in rows:
            for source in str(row["natural_product_sources"]).split(";"):
                if source:
                    latest_source_counts[source] += 1
        eligible = len(rows) >= MIN_LATEST_MOLECULES and natural_count >= MIN_LATEST_NATURAL_PRODUCTS
        target_result: dict[str, Any] = {
            "target_chembl_id": target_id,
            "target_name": spec["target_name"],
            "status": "included" if eligible else "excluded_insufficient_temporal_or_natural_product_coverage",
            "selected_context_key": context_key,
            "selected_context": selected[0].get("assay_context") or {},
            "normalized_measurements": len(all_records),
            "compatible_measurements": len(selected),
            "compatible_molecules": len({row["molecule_chembl_id"] for row in selected}),
            "latest_publication_molecules": len(rows),
            "latest_publication_natural_products": natural_count,
            "latest_natural_product_source_molecules": dict(sorted(latest_source_counts.items())),
            "exclusions": {"context": context_exclusions, "temporal": split_exclusions},
        }
        input_manifest.append({
            "target_chembl_id": target_id,
            "activities": str(activities_path.resolve()),
            "activities_sha256": _sha256(activities_path),
            "natural_products": str(natural_path.resolve()),
            "natural_products_sha256": _sha256(natural_path),
        })
        if not eligible:
            targets[target_id] = target_result
            continue

        development, evaluation = deterministic_partition(
            rows, development_fraction=development_fraction, seed=seed + target_index
        )
        partition_rows[target_id] = {"development": development, "evaluation": evaluation}
        for partition_name, current_rows in (("development", development), ("evaluation", evaluation)):
            for row in current_rows:
                all_cohort_rows.append({
                    "target_chembl_id": target_id,
                    "target_name": spec["target_name"],
                    "partition": partition_name,
                    **row,
                })
        target_result["partitions"] = {}
        for partition_name, current_rows in (("development", development), ("evaluation", evaluation)):
            target_result["partitions"][partition_name] = {
                "all": _metric_table(current_rows, k=k, active_threshold=active_threshold),
                "natural_products": _metric_table(
                    [row for row in current_rows if row["is_verified_natural_product"]],
                    k=k,
                    active_threshold=active_threshold,
                ),
            }
            component_rows = component_predictive_value(
                current_rows, iterations=iterations, seed=seed + target_index * 100 + (0 if partition_name == "development" else 50)
            )
            for component in component_rows:
                all_component_rows.append({
                    "target_chembl_id": target_id,
                    "partition": partition_name,
                    "component": component["component"],
                    "n_molecules": component["n_molecules"],
                    "signed_residual_rho": component["signed_residual"]["spearman_rho"],
                    "signed_residual_ci_low": component["signed_residual"]["ci95"][0],
                    "signed_residual_ci_high": component["signed_residual"]["ci95"][1],
                    "signed_residual_p": component["signed_residual"]["permutation_p"],
                    "absolute_error_rho": component["absolute_error"]["spearman_rho"],
                    "absolute_error_ci_low": component["absolute_error"]["ci95"][0],
                    "absolute_error_ci_high": component["absolute_error"]["ci95"][1],
                    "absolute_error_p": component["absolute_error"]["permutation_p"],
                })

        error_rows, harmful_pairs, error_summary = compound_error_analysis(evaluation)
        target_result["evaluation_error_analysis"] = error_summary
        for row in error_rows:
            all_error_rows.append({"target_chembl_id": target_id, **row})
        for row in harmful_pairs:
            all_harmful_pairs.append({"target_chembl_id": target_id, **row})

        development_ids = {str(row["molecule_id"]) for row in development}
        evaluation_ids = {str(row["molecule_id"]) for row in evaluation}
        for partition_name, ids in (("development", development_ids), ("evaluation", evaluation_ids)):
            for row in assay_transferability(all_records, ids):
                all_transfer_rows.append({"target_chembl_id": target_id, "partition": partition_name, **row})
        targets[target_id] = target_result
        dataset_rows.append({
            "target_chembl_id": target_id,
            "target_name": spec["target_name"],
            "status": target_result["status"],
            "normalized_measurements": target_result["normalized_measurements"],
            "compatible_measurements": target_result["compatible_measurements"],
            "compatible_molecules": target_result["compatible_molecules"],
            "latest_publication_molecules": target_result["latest_publication_molecules"],
            "latest_publication_natural_products": target_result["latest_publication_natural_products"],
            "latest_npass_molecules": latest_source_counts.get("NPASS", 0),
            "latest_coconut_molecules": latest_source_counts.get("COCONUT", 0),
        })

    included_ids = [target_id for target_id, target in targets.items() if target["status"] == "included"]
    if len(included_ids) < 3:
        raise ValueError(f"Only {len(included_ids)} targets met prespecified coverage; at least 3 are required")
    development_rows = {target_id: partition_rows[target_id]["development"] for target_id in included_ids}
    evaluation_rows = {target_id: partition_rows[target_id]["evaluation"] for target_id in included_ids}
    development_comparisons = {
        method: hierarchical_macro_inference(
            development_rows,
            method,
            k=k,
            active_threshold=active_threshold,
            iterations=iterations,
            seed=seed + 10_000 + method_index,
        )
        for method_index, method in enumerate(PRIMARY_ALTERNATIVES)
    }
    selection = select_on_development(development_comparisons)
    evaluation_comparisons = {
        method: hierarchical_macro_inference(
            evaluation_rows,
            method,
            k=k,
            active_threshold=active_threshold,
            iterations=iterations,
            seed=seed + 20_000 + method_index,
        )
        for method_index, method in enumerate(PRIMARY_ALTERNATIVES)
    }
    selected_method = selection["selected_method"]
    evaluation_supports_superiority = False
    if selected_method != "median_potency":
        selected_evaluation = evaluation_comparisons[selected_method]
        evaluation_supports_superiority = all(
            float(selected_evaluation[metric]["difference_ci95"][0]) > 0
            and float(selected_evaluation[metric]["paired_randomization_p"]) < 0.05
            for metric in ("ndcg_at_k", "kendall_tau")
        )

    connection = sqlite3.connect(natural_product_database)
    index_summary = {
        "unique_structures": connection.execute("SELECT COUNT(*) FROM compounds").fetchone()[0],
        "provenance_records": connection.execute("SELECT COUNT(*) FROM provenance").fetchone()[0],
        "npass_records": connection.execute("SELECT COUNT(*) FROM provenance WHERE source='NPASS'").fetchone()[0],
        "coconut_records": connection.execute("SELECT COUNT(*) FROM provenance WHERE source='COCONUT'").fetchone()[0],
        "rejected_records": connection.execute("SELECT COUNT(*) FROM rejects").fetchone()[0],
        "structures_with_both_npass_and_coconut": connection.execute(
            """SELECT COUNT(*) FROM (
                   SELECT compound_id FROM provenance
                   WHERE source IN ('NPASS', 'COCONUT')
                   GROUP BY compound_id HAVING COUNT(DISTINCT source) = 2
               )"""
        ).fetchone()[0],
        "source_versions": {
            row[0]: sorted(value for value in row[1].split(",") if value)
            for row in connection.execute(
                "SELECT source, GROUP_CONCAT(DISTINCT dataset_version) FROM provenance GROUP BY source"
            )
        },
        "sha256": _sha256(natural_product_database),
    }
    connection.close()
    results: dict[str, Any] = {
        "experiment": "PhytoTruth Multi-Target Evidence Benchmark",
        "generated_at": datetime.now(UTC).isoformat(),
        "selection_was_frozen_before_evaluation": True,
        "seed": seed,
        "iterations": iterations,
        "k": k,
        "active_threshold": active_threshold,
        "development_fraction": development_fraction,
        "eligibility": {
            "minimum_latest_publication_molecules": MIN_LATEST_MOLECULES,
            "minimum_latest_publication_natural_products": MIN_LATEST_NATURAL_PRODUCTS,
            "context_selection": "maximum distinct molecule coverage; no activity-value criterion",
        },
        "natural_product_index": index_summary,
        "targets": targets,
        "included_targets": included_ids,
        "development_macro_comparisons": development_comparisons,
        "development_selection": selection,
        "evaluation_macro_comparisons": evaluation_comparisons,
        "scientific_conclusion": {
            "selected_ranking": selected_method,
            "superiority_claim": evaluation_supports_superiority,
            "status": (
                "independent_superiority_supported"
                if evaluation_supports_superiority
                else "superiority_not_demonstrated"
            ),
            "organism_transfer": "not_estimable_all_selected_targets_and_available_contexts_are_human",
        },
    }
    _write_csv(output_dir / "dataset_summary.csv", dataset_rows)
    _write_csv(output_dir / "benchmark_cohort.csv", all_cohort_rows)
    _write_csv(output_dir / "metric_summary.csv", _flatten_metrics(results))
    _write_csv(output_dir / "component_predictive_value.csv", all_component_rows)
    _write_csv(output_dir / "compound_error_analysis.csv", all_error_rows)
    _write_csv(output_dir / "harmful_reversals.csv", all_harmful_pairs)
    _write_csv(output_dir / "assay_transferability.csv", all_transfer_rows)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    _comparison_svg(output_dir / "comparison.svg", results)
    manifest = {
        "specification": str(specs_path.resolve()),
        "specification_sha256": _sha256(specs_path),
        "natural_product_database": str(natural_product_database.resolve()),
        "natural_product_database_sha256": _sha256(natural_product_database),
        "inputs": input_manifest,
        "results_sha256": _sha256(output_dir / "results.json"),
        "leakage_contract": "temporal evidence split precedes molecule development/evaluation split; evaluation is not used for selection",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    command = (
        "python -m app.multi_target_benchmark "
        f'--specs "{specs_path}" --natural-product-database "{natural_product_database}" '
        f'--output-dir "{output_dir}" --development-fraction {development_fraction} --k {k} '
        f"--active-threshold {active_threshold} --iterations {iterations} --seed {seed}"
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the leakage-resistant multi-target evidence benchmark")
    parser.add_argument("--specs", type=Path, required=True)
    parser.add_argument("--natural-product-database", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-fraction", type=float, default=0.6)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--active-threshold", type=float, default=6.0)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    results = run_benchmark(
        args.specs,
        args.natural_product_database,
        args.output_dir,
        development_fraction=args.development_fraction,
        k=args.k,
        active_threshold=args.active_threshold,
        iterations=args.iterations,
        seed=args.seed,
    )
    print(json.dumps({
        "included_targets": results["included_targets"],
        "selected_ranking": results["development_selection"]["selected_method"],
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
