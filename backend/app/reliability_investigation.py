from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import UTC, datetime
from html import escape
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from scipy.stats import rankdata

from .reversal import METRIC_NAMES, evaluate_ranking, paired_ranking_inference
from .reversal_experiment import (
    _jsonl,
    _sha256,
    build_benchmark_rows,
    latest_publication_split,
    natural_product_membership,
    select_compatible_measurements,
)


STRONG_BASELINES = (
    "baseline_max_potency",
    "median_potency",
    "mean_potency",
    "latest_training_potency",
    "lower_confidence_bound",
)
ABLATIONS = (
    "ablate_reliability",
    "ablate_uncertainty",
    "ablate_explicit_consistency",
    "ablate_support",
    "phytotruth_full",
)
ALL_METHODS = STRONG_BASELINES + ABLATIONS


def deterministic_partition(
    rows: Sequence[dict[str, Any]],
    *,
    development_fraction: float = 0.6,
    seed: int = 2026,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Stratify by NP membership using only molecule IDs, never outcomes or scores."""
    if not 0 < development_fraction < 1:
        raise ValueError("development_fraction must be between zero and one")
    molecule_ids = [str(row["molecule_id"]) for row in rows]
    if len(molecule_ids) != len(set(molecule_ids)):
        raise ValueError("molecule IDs must be unique before partitioning")

    groups: dict[bool, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[bool(row["is_verified_natural_product"])].append(dict(row))

    development: list[dict[str, Any]] = []
    evaluation: list[dict[str, Any]] = []
    for is_natural in sorted(groups):
        group = sorted(
            groups[is_natural],
            key=lambda row: hashlib.sha256(
                f"{seed}:{row['molecule_id']}".encode("utf-8")
            ).hexdigest(),
        )
        if len(group) == 1:
            evaluation.extend(group)
            continue
        n_development = min(len(group) - 1, max(1, round(len(group) * development_fraction)))
        development.extend(group[:n_development])
        evaluation.extend(group[n_development:])
    return sorted(development, key=lambda row: row["molecule_id"]), sorted(
        evaluation, key=lambda row: row["molecule_id"]
    )


def _metric_table(rows: Sequence[dict[str, Any]], *, k: int, active_threshold: float) -> dict[str, Any]:
    truth = [float(row["heldout_median_pactivity"]) for row in rows]
    actual_k = min(k, len(rows))
    return {
        "n_molecules": len(rows),
        "verified_natural_products": sum(bool(row["is_verified_natural_product"]) for row in rows),
        "heldout_active_molecules": sum(value >= active_threshold for value in truth),
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


def _inference(
    rows: Sequence[dict[str, Any]],
    candidate: str,
    *,
    k: int,
    active_threshold: float,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    if len(rows) < 2:
        return {"status": "insufficient_molecules"}
    return paired_ranking_inference(
        [float(row["median_potency"]) for row in rows],
        [float(row[candidate]) for row in rows],
        [float(row["heldout_median_pactivity"]) for row in rows],
        k=min(k, len(rows)),
        active_threshold=active_threshold,
        bootstrap_iterations=iterations,
        permutation_iterations=iterations,
        seed=seed,
    )


def select_ranking_on_development(inference: dict[str, Any]) -> dict[str, Any]:
    """Use a prespecified conservative gate; evaluation outcomes are not an input."""
    supported = all(
        float(inference[metric]["difference_ci95"][0]) > 0
        and float(inference[metric]["paired_randomization_p"]) < 0.05
        for metric in ("ndcg_at_k", "kendall_tau")
    )
    return {
        "selected_method": "phytotruth_full" if supported else "median_potency",
        "gate": (
            "Select the current evidence adjustment only when development-set nDCG@K and "
            "Kendall tau-b both improve over median potency, both paired 95% interval lower "
            "bounds exceed zero, and both two-sided score-swap p-values are below 0.05."
        ),
        "passed": supported,
    }


def _spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3 or np.ptp(left) == 0 or np.ptp(right) == 0:
        return 0.0
    return float(np.corrcoef(rankdata(left), rankdata(right))[0, 1])


def _correlation_inference(
    signal: Sequence[float],
    outcome: Sequence[float],
    *,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    x = np.asarray(signal, dtype=float)
    y = np.asarray(outcome, dtype=float)
    point = _spearman(x, y)
    rng = np.random.default_rng(seed)
    bootstrap: list[float] = []
    for _ in range(iterations):
        indices = rng.integers(0, len(x), size=len(x))
        bootstrap.append(_spearman(x[indices], y[indices]))
    exceedances = 0
    for _ in range(iterations):
        permuted = rng.permutation(y)
        exceedances += int(abs(_spearman(x, permuted)) >= abs(point) - 1e-12)
    return {
        "spearman_rho": point,
        "ci95": [float(value) for value in np.quantile(bootstrap, [0.025, 0.975])],
        "permutation_p": (exceedances + 1) / (iterations + 1),
    }


def component_predictive_value(
    rows: Sequence[dict[str, Any]], *, iterations: int = 2000, seed: int = 2026
) -> list[dict[str, Any]]:
    """Test whether pre-holdout components predict signed drift or future absolute error."""
    specifications = {
        "assay_support": (
            lambda row: math.log1p(float(row["train_assays"])),
            lambda row: -math.log1p(float(row["train_assays"])),
        ),
        "independent_study_support": (
            lambda row: math.log1p(float(row["train_documents"])),
            lambda row: -math.log1p(float(row["train_documents"])),
        ),
        "experimental_consistency": (
            lambda row: -float(row["consistency_penalty"]),
            lambda row: float(row["consistency_penalty"]),
        ),
        "bootstrap_uncertainty": (
            lambda row: -float(row["uncertainty_penalty"]),
            lambda row: float(row["uncertainty_penalty"]),
        ),
        "composite_reliability": (
            lambda row: -float(row["reliability_penalty"]),
            lambda row: float(row["reliability_penalty"]),
        ),
        "current_total_adjustment": (
            lambda row: float(row["phytotruth_full"]) - float(row["median_potency"]),
            lambda row: abs(float(row["phytotruth_full"]) - float(row["median_potency"])),
        ),
    }
    residual = [float(row["heldout_median_pactivity"]) - float(row["median_potency"]) for row in rows]
    absolute_error = [abs(value) for value in residual]
    results: list[dict[str, Any]] = []
    for index, (name, (adjustment, risk)) in enumerate(specifications.items()):
        signed = _correlation_inference(
            [adjustment(row) for row in rows], residual, iterations=iterations, seed=seed + index * 2
        )
        error = _correlation_inference(
            [risk(row) for row in rows], absolute_error, iterations=iterations, seed=seed + index * 2 + 1
        )
        results.append({
            "component": name,
            "n_molecules": len(rows),
            "signed_residual": signed,
            "absolute_error": error,
        })
    return results


def _component_adjustments(row: dict[str, Any]) -> dict[str, float]:
    return {
        "assay_support": float(row["support_bonus"]),
        "experimental_consistency": -float(row["consistency_penalty"]),
        "composite_reliability": -float(row["reliability_penalty"]),
        "bootstrap_uncertainty": -float(row["uncertainty_penalty"]),
    }


def compound_error_analysis(rows: Sequence[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Identify pairwise order changes where the heuristic breaks a correct median order."""
    rank_fields: dict[str, dict[str, int]] = {}
    for method in ("median_potency", "phytotruth_full", "heldout_median_pactivity"):
        ordered = sorted(rows, key=lambda row: (-float(row[method]), str(row["molecule_id"])))
        rank_fields[method] = {str(row["molecule_id"]): index + 1 for index, row in enumerate(ordered)}

    compounds: dict[str, dict[str, Any]] = {}
    for row in rows:
        molecule_id = str(row["molecule_id"])
        compounds[molecule_id] = {
            "molecule_id": molecule_id,
            "is_verified_natural_product": bool(row["is_verified_natural_product"]),
            "median_rank": rank_fields["median_potency"][molecule_id],
            "phytotruth_rank": rank_fields["phytotruth_full"][molecule_id],
            "heldout_rank": rank_fields["heldout_median_pactivity"][molecule_id],
            "median_potency": float(row["median_potency"]),
            "phytotruth_full": float(row["phytotruth_full"]),
            "heldout_median_pactivity": float(row["heldout_median_pactivity"]),
            "score_adjustment": float(row["phytotruth_full"]) - float(row["median_potency"]),
            "harmful_pair_count": 0,
            "beneficial_pair_count": 0,
        }

    harmful_pairs: list[dict[str, Any]] = []
    dominant_harms: Counter[str] = Counter()
    beneficial = 0
    comparable = 0
    for first_index, first in enumerate(rows):
        for second in rows[first_index + 1:]:
            truth_delta = float(first["heldout_median_pactivity"]) - float(second["heldout_median_pactivity"])
            median_delta = float(first["median_potency"]) - float(second["median_potency"])
            full_delta = float(first["phytotruth_full"]) - float(second["phytotruth_full"])
            if truth_delta == 0 or median_delta == 0 or full_delta == 0:
                continue
            comparable += 1
            median_correct = math.copysign(1, truth_delta) == math.copysign(1, median_delta)
            full_correct = math.copysign(1, truth_delta) == math.copysign(1, full_delta)
            if median_correct == full_correct:
                continue
            first_id = str(first["molecule_id"])
            second_id = str(second["molecule_id"])
            if not median_correct and full_correct:
                beneficial += 1
                compounds[first_id]["beneficial_pair_count"] += 1
                compounds[second_id]["beneficial_pair_count"] += 1
                continue
            deltas = {
                name: _component_adjustments(first)[name] - _component_adjustments(second)[name]
                for name in _component_adjustments(first)
            }
            opposing = {
                name: value
                for name, value in deltas.items()
                if value != 0 and math.copysign(1, value) != math.copysign(1, median_delta)
            }
            dominant = max(opposing, key=lambda name: abs(opposing[name])) if opposing else "combined_or_tie"
            dominant_harms[dominant] += 1
            compounds[first_id]["harmful_pair_count"] += 1
            compounds[second_id]["harmful_pair_count"] += 1
            harmful_pairs.append({
                "molecule_a": first_id,
                "molecule_b": second_id,
                "natural_product_pair": bool(first["is_verified_natural_product"] and second["is_verified_natural_product"]),
                "truth_delta": truth_delta,
                "median_delta": median_delta,
                "phytotruth_delta": full_delta,
                "dominant_harmful_component": dominant,
                **{f"{name}_pair_adjustment": value for name, value in deltas.items()},
            })
    compound_rows = sorted(
        compounds.values(),
        key=lambda row: (-int(row["harmful_pair_count"]), str(row["molecule_id"])),
    )
    harmful_pairs.sort(key=lambda row: (-abs(float(row["truth_delta"])), row["molecule_a"], row["molecule_b"]))
    summary = {
        "comparable_pairs": comparable,
        "harmful_reversals": len(harmful_pairs),
        "beneficial_reversals": beneficial,
        "net_harmful_reversals": len(harmful_pairs) - beneficial,
        "dominant_component_counts": dict(sorted(dominant_harms.items())),
    }
    return compound_rows, harmful_pairs, summary


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _flatten_component_rows(partition: str, subgroup: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    flattened: list[dict[str, Any]] = []
    for row in rows:
        flattened.append({
            "partition": partition,
            "subgroup": subgroup,
            "component": row["component"],
            "n_molecules": row["n_molecules"],
            "signed_residual_rho": row["signed_residual"]["spearman_rho"],
            "signed_residual_ci_low": row["signed_residual"]["ci95"][0],
            "signed_residual_ci_high": row["signed_residual"]["ci95"][1],
            "signed_residual_p": row["signed_residual"]["permutation_p"],
            "absolute_error_rho": row["absolute_error"]["spearman_rho"],
            "absolute_error_ci_low": row["absolute_error"]["ci95"][0],
            "absolute_error_ci_high": row["absolute_error"]["ci95"][1],
            "absolute_error_p": row["absolute_error"]["permutation_p"],
        })
    return flattened


def _metric_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for partition in ("development", "evaluation"):
        for subgroup in ("all", "natural_products"):
            table = results["partitions"][partition][subgroup]["metric_table"]
            for method, metrics in table["metrics"].items():
                for metric, value in metrics.items():
                    output.append({
                        "partition": partition,
                        "subgroup": subgroup,
                        "method": method,
                        "metric": metric,
                        "value": value,
                        "n_molecules": table["n_molecules"],
                        "k": table["k"],
                    })
    return output


def _comparison_svg(path: Path, results: dict[str, Any]) -> None:
    metrics = results["partitions"]["evaluation"]["all"]["metric_table"]["metrics"]
    methods = ("baseline_max_potency", "median_potency", "phytotruth_full")
    labels = ("Maximum potency", "Median potency", "Current PhytoTruth")
    colours = ("#98a2b3", "#12b76a", "#f04438")
    width, height = 820, 390
    chart_left, chart_top, chart_width, chart_height = 90, 55, 650, 250
    bars: list[str] = []
    for metric_index, metric in enumerate(("ndcg_at_k", "kendall_tau")):
        group_x = chart_left + metric_index * 330
        for index, method in enumerate(methods):
            value = float(metrics[method][metric])
            bar_height = max(0.0, value) * chart_height
            x = group_x + index * 62
            y = chart_top + chart_height - bar_height
            bars.append(f'<rect x="{x}" y="{y:.1f}" width="48" height="{bar_height:.1f}" rx="5" fill="{colours[index]}"/>')
            bars.append(f'<text x="{x + 24}" y="{y - 7:.1f}" text-anchor="middle" font-size="12" fill="#344054">{value:.3f}</text>')
        bars.append(f'<text x="{group_x + 86}" y="{chart_top + chart_height + 28}" text-anchor="middle" font-size="14" fill="#101828">{escape(metric.replace("_", " "))}</text>')
    legend = "".join(
        f'<rect x="{150 + index * 190}" y="350" width="12" height="12" rx="2" fill="{colours[index]}"/><text x="{168 + index * 190}" y="361" font-size="12" fill="#475467">{escape(label)}</text>'
        for index, label in enumerate(labels)
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#ffffff"/><text x="36" y="30" font-size="18" font-family="sans-serif" fill="#101828">Untouched evaluation partition</text>
<line x1="{chart_left - 10}" y1="{chart_top + chart_height}" x2="{chart_left + chart_width}" y2="{chart_top + chart_height}" stroke="#d0d5dd"/>
<g font-family="sans-serif">{''.join(bars)}{legend}</g></svg>'''
    path.write_text(svg, encoding="utf-8")


def _report(results: dict[str, Any]) -> str:
    dev = results["partitions"]["development"]["all"]
    evaluation = results["partitions"]["evaluation"]["all"]
    natural = results["partitions"]["evaluation"]["natural_products"]
    selection = results["development_selection"]
    eval_metrics = evaluation["metric_table"]["metrics"]
    natural_metrics = natural["metric_table"]["metrics"]
    current_inference = evaluation["comparisons_to_median"]["phytotruth_full"]
    natural_inference = natural["comparisons_to_median"]["phytotruth_full"]
    pairs = results["error_analysis"]["evaluation_all"]
    component_counts = ", ".join(
        f"{key}: {value}" for key, value in pairs["dominant_component_counts"].items()
    ) or "none"
    return f"""# Evidence Reliability Engine scientific investigation

Generated {results['generated_at']} from normalized ChEMBL measurements and the recorded NPASS structure cross-match.

## Prespecified validation design

The latest-publication-per-molecule protocol was retained: all score inputs precede each molecule's held-out publication evidence, and documents are disjoint across time. Molecules were then deterministically partitioned, using only a hash of molecule ID and NPASS membership, into {dev['metric_table']['n_molecules']} development and {evaluation['metric_table']['n_molecules']} untouched-for-revision evaluation molecules. The evaluation outcomes were not inputs to the selection gate. This partition is untouched for this revision procedure, but it is not an external validation dataset because aggregate results on the parent cohort were reported previously.

## Why the current heuristic fails

The current score subtracts inconsistency once in `evidence_score`, then adds a composite reliability penalty that includes the same MAD/contradiction signal, and finally subtracts bootstrap interval width. This double/triple use of dispersion is not calibrated to future assay outcomes. Its positive support term counts assays, while independent publications are only weakly and saturatingly represented inside the composite reliability value. Consequently, the heuristic can reverse a robust median ordering when sparse or heterogeneous evidence produces a large penalty.

On the untouched evaluation partition, current PhytoTruth minus median potency Kendall tau-b was {current_inference['kendall_tau']['difference']:.3f} (95% CI {current_inference['kendall_tau']['difference_ci95'][0]:.3f} to {current_inference['kendall_tau']['difference_ci95'][1]:.3f}; p={current_inference['kendall_tau']['paired_randomization_p']:.4f}). nDCG@10 difference was {current_inference['ndcg_at_k']['difference']:.3f} (95% CI {current_inference['ndcg_at_k']['difference_ci95'][0]:.3f} to {current_inference['ndcg_at_k']['difference_ci95'][1]:.3f}).

Pairwise error analysis found {pairs['harmful_reversals']} harmful and {pairs['beneficial_reversals']} beneficial changes to median-correct orderings among {pairs['comparable_pairs']} comparable evaluation pairs. Dominant harmful components were {component_counts}. See `compound_error_analysis.csv` and `harmful_reversals.csv` for traceable molecule-level details.

## Natural products

The untouched NPASS-verified evaluation subgroup has only {natural['metric_table']['n_molecules']} molecules. Current PhytoTruth minus median nDCG@{natural['metric_table']['k']} was {natural_inference['ndcg_at_k']['difference']:.3f} (95% CI {natural_inference['ndcg_at_k']['difference_ci95'][0]:.3f} to {natural_inference['ndcg_at_k']['difference_ci95'][1]:.3f}); Kendall difference was {natural_inference['kendall_tau']['difference']:.3f}. Median nDCG was {natural_metrics['median_potency']['ndcg_at_k']:.3f}, versus {natural_metrics['phytotruth_full']['ndcg_at_k']:.3f} for the current heuristic. This sample is too small for a reliable subgroup superiority claim. COCONUT did not contribute verified compounds to this target artifact; this subgroup is NPASS-only.

## Development-only decision and defensible change

The prespecified gate was: {selection['gate']} It did {'pass' if selection['passed'] else 'not pass'} on development data, so the selected ranking is `{selection['selected_method']}`. No new weights were fitted. The scientifically defensible change is therefore a conservative guardrail: use robust median potency for prioritization and retain reliability fields as descriptive warnings, not as validated ranking probabilities or penalties.

Untouched evaluation median metrics were nDCG@10={eval_metrics['median_potency']['ndcg_at_k']:.3f}, precision@10={eval_metrics['median_potency']['precision_at_k']:.3f}, enrichment={eval_metrics['median_potency']['enrichment_factor_at_k']:.3f}, and Kendall tau-b={eval_metrics['median_potency']['kendall_tau']:.3f}. Current-heuristic metrics were nDCG@10={eval_metrics['phytotruth_full']['ndcg_at_k']:.3f}, precision@10={eval_metrics['phytotruth_full']['precision_at_k']:.3f}, enrichment={eval_metrics['phytotruth_full']['enrichment_factor_at_k']:.3f}, and Kendall tau-b={eval_metrics['phytotruth_full']['kendall_tau']:.3f}.

## Predictive-value and ablation interpretation

`component_predictive_value.csv` reports development and evaluation Spearman associations, molecule-bootstrap 95% intervals, and permutation p-values for two questions: whether a component points in the direction of later potency drift, and whether its risk signal predicts absolute future error. These are association tests, not probability calibration. `metric_summary.csv` compares maximum, mean, most-recent, lower-confidence-bound, and median potency with every existing ablation.

## Limitations and required data

- The compatible ChEMBL context is one target/context; generalization to other targets is unknown.
- The revision evaluation partition contains {evaluation['metric_table']['n_molecules']} molecules and only {natural['metric_table']['n_molecules']} NPASS-verified natural products, so natural-product confidence intervals are wide.
- Publication identifiers are a proxy for independent studies; laboratory, protocol, batch, and replicate identifiers are unavailable.
- Bootstrap intervals over repeated measurements do not capture systematic assay bias.
- A scientifically credible learned reliability adjustment needs a larger multi-target dataset, externally held-out publications/laboratories, richer assay-protocol covariates, and substantially more COCONUT/NPASS-verified molecules. Until then, reliability indicators should remain explanatory warnings.
"""


def run_investigation(
    activities: Path,
    natural_products: Path,
    output_dir: Path,
    *,
    context_key: str | None = None,
    development_fraction: float = 0.6,
    k: int = 10,
    active_threshold: float = 6.0,
    iterations: int = 2000,
    seed: int = 2026,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    natural_members = natural_product_membership(natural_products)
    selected, selected_context, context_exclusions = select_compatible_measurements(
        _jsonl(activities), context_key=context_key
    )
    latest, split_exclusions = latest_publication_split(selected)
    rows = build_benchmark_rows(latest, natural_members)
    development, evaluation = deterministic_partition(
        rows, development_fraction=development_fraction, seed=seed
    )

    def analyse_partition(partition_rows: list[dict[str, Any]], partition_seed: int) -> dict[str, Any]:
        subgroups = {
            "all": partition_rows,
            "natural_products": [row for row in partition_rows if row["is_verified_natural_product"]],
        }
        output: dict[str, Any] = {}
        for subgroup_index, (name, subgroup) in enumerate(subgroups.items()):
            output[name] = {
                "metric_table": _metric_table(subgroup, k=k, active_threshold=active_threshold),
                "comparisons_to_median": {
                    method: _inference(
                        subgroup,
                        method,
                        k=k,
                        active_threshold=active_threshold,
                        iterations=iterations,
                        seed=partition_seed + subgroup_index * 100 + method_index,
                    )
                    for method_index, method in enumerate(ALL_METHODS)
                    if method != "median_potency"
                },
                "component_predictive_value": component_predictive_value(
                    subgroup, iterations=iterations, seed=partition_seed + subgroup_index * 1000
                ) if len(subgroup) >= 3 else [],
            }
        return output

    development_results = analyse_partition(development, seed + 10_000)
    development_selection = select_ranking_on_development(
        development_results["all"]["comparisons_to_median"]["phytotruth_full"]
    )
    evaluation_results = analyse_partition(evaluation, seed + 20_000)

    development_compounds, development_pairs, development_pair_summary = compound_error_analysis(development)
    evaluation_compounds, evaluation_pairs, evaluation_pair_summary = compound_error_analysis(evaluation)
    evaluation_natural_rows = [row for row in evaluation if row["is_verified_natural_product"]]
    natural_compounds, natural_pairs, natural_pair_summary = compound_error_analysis(evaluation_natural_rows)

    context = selected[0].get("assay_context") or {}
    results: dict[str, Any] = {
        "experiment": "PhytoTruth Evidence Reliability Investigation",
        "generated_at": datetime.now(UTC).isoformat(),
        "target_chembl_id": context.get("target_chembl_id"),
        "assay_context_key": selected_context,
        "assay_context": context,
        "selected_measurements": len(selected),
        "selected_molecules": len({record["molecule_chembl_id"] for record in selected}),
        "latest_publication_cohort_molecules": len(rows),
        "latest_publication_np_molecules": sum(bool(row["is_verified_natural_product"]) for row in rows),
        "natural_product_verification_sources": sorted({
            source
            for row in rows
            for source in str(row["natural_product_sources"]).split(";")
            if source
        }),
        "development_fraction": development_fraction,
        "seed": seed,
        "k": k,
        "active_threshold": active_threshold,
        "iterations": iterations,
        "partition_contract": "stratified deterministic SHA-256 of seed:molecule_id; no assay outcomes or scores",
        "temporal_contract": "latest publication year held out; training publications are strictly earlier and document-disjoint",
        "exclusions": {"context": context_exclusions, "latest_publication": split_exclusions},
        "development_selection": development_selection,
        "partitions": {"development": development_results, "evaluation": evaluation_results},
        "error_analysis": {
            "development_all": development_pair_summary,
            "evaluation_all": evaluation_pair_summary,
            "evaluation_natural_products": natural_pair_summary,
        },
        "scientific_conclusion": {
            "status": "current_adjustment_not_validated" if not development_selection["passed"] else "development_gate_passed",
            "selected_ranking": development_selection["selected_method"],
            "reliability_use": "descriptive_warning_only" if not development_selection["passed"] else "ranking_candidate_requires_external_replication",
        },
    }

    for row in development:
        row["revision_partition"] = "development"
    for row in evaluation:
        row["revision_partition"] = "evaluation"
    _write_csv(output_dir / "development_cohort.csv", development)
    _write_csv(output_dir / "evaluation_cohort.csv", evaluation)
    _write_csv(output_dir / "compound_error_analysis_development.csv", development_compounds)
    _write_csv(output_dir / "compound_error_analysis.csv", evaluation_compounds)
    _write_csv(output_dir / "compound_error_analysis_natural_products.csv", natural_compounds)
    _write_csv(output_dir / "harmful_reversals_development.csv", development_pairs)
    _write_csv(output_dir / "harmful_reversals.csv", evaluation_pairs)
    _write_csv(output_dir / "harmful_reversals_natural_products.csv", natural_pairs)
    component_rows: list[dict[str, Any]] = []
    for partition in ("development", "evaluation"):
        for subgroup in ("all", "natural_products"):
            component_rows.extend(_flatten_component_rows(
                partition,
                subgroup,
                results["partitions"][partition][subgroup]["component_predictive_value"],
            ))
    _write_csv(output_dir / "component_predictive_value.csv", component_rows)
    _write_csv(output_dir / "metric_summary.csv", _metric_rows(results))
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    _comparison_svg(output_dir / "comparison.svg", results)
    command = (
        "python -m app.reliability_investigation "
        f'--activities "{activities}" --natural-products "{natural_products}" '
        f'--output-dir "{output_dir}" --development-fraction {development_fraction} '
        f"--k {k} --active-threshold {active_threshold} --iterations {iterations} --seed {seed}"
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    manifest = {
        "activities": str(activities.resolve()),
        "activities_sha256": _sha256(activities),
        "natural_products": str(natural_products.resolve()),
        "natural_products_sha256": _sha256(natural_products),
        "results_sha256": _sha256(output_dir / "results.json"),
        "partition_sha256": hashlib.sha256(
            "\n".join(
                f"development:{row['molecule_id']}" for row in development
            ).encode("utf-8") + "\n".join(
                f"evaluation:{row['molecule_id']}" for row in evaluation
            ).encode("utf-8")
        ).hexdigest(),
        "code_contract": "evaluation outcomes are not inputs to development selection",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Investigate the Evidence Reliability Engine without test leakage")
    parser.add_argument("--activities", type=Path, required=True)
    parser.add_argument("--natural-products", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-key")
    parser.add_argument("--development-fraction", type=float, default=0.6)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--active-threshold", type=float, default=6.0)
    parser.add_argument("--iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    results = run_investigation(
        args.activities,
        args.natural_products,
        args.output_dir,
        context_key=args.context_key,
        development_fraction=args.development_fraction,
        k=args.k,
        active_threshold=args.active_threshold,
        iterations=args.iterations,
        seed=args.seed,
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "development_molecules": results["partitions"]["development"]["all"]["metric_table"]["n_molecules"],
        "evaluation_molecules": results["partitions"]["evaluation"]["all"]["metric_table"]["n_molecules"],
        "selected_ranking": results["development_selection"]["selected_method"],
    }, indent=2))


if __name__ == "__main__":
    main()
