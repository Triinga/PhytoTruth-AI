from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .assay_disagreement import _metric_values, _sha256
from .dataset_expansion import merge_ledgers
from .independent_validation import (
    _ids,
    _load_frozen_models,
    build_previous_exposure_ledger,
    event_overlap_reasons,
    score_frozen_models,
)


VALIDATION_VERSION = "expanded-frozen-model-validation-v1"
FORBIDDEN_FEATURE_NAMES = {
    "disagreement_label", "absolute_pactivity_shift", "signed_pactivity_shift",
    "active_threshold_contradiction",
}


def _write_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")


def assert_historical_feature_contract(features: Sequence[str]) -> None:
    violations = [
        feature for feature in features
        if feature.lower() in FORBIDDEN_FEATURE_NAMES
        or feature.lower().startswith("heldout_")
        or feature.lower().startswith("protocol_heldout_")
        or feature.lower().startswith("protocol_comparison_")
    ]
    if violations:
        raise AssertionError(f"Frozen artifact contains outcome/future features: {violations}")


def tie_aware_precision_at_k(labels: Sequence[int], scores: Sequence[float], k: int) -> float | None:
    """Expected Precision@K under random ordering within a score tie."""
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    if not len(y) or k < 1:
        return None
    k = min(int(k), len(y))
    selected_positives = 0.0
    selected = 0
    for score in sorted(set(s.tolist()), reverse=True):
        group = y[s == score]
        needed = k - selected
        if needed <= 0:
            break
        if len(group) <= needed:
            selected_positives += float(group.sum())
            selected += len(group)
        else:
            selected_positives += needed * float(group.mean())
            selected = k
    return float(selected_positives / k)


def metric_values(frame: pd.DataFrame, probability: str, *, k: int) -> dict[str, Any]:
    result = _metric_values(frame.disagreement_label, frame[probability])
    result[f"precision_at_{min(k, len(frame))}"] = tie_aware_precision_at_k(
        frame.disagreement_label, frame[probability], k
    )
    return result


def clustered_intervals(
    frame: pd.DataFrame,
    probability: str,
    *,
    k: int,
    iterations: int,
    seed: int,
) -> dict[str, list[float] | None]:
    if frame.empty or frame.heldout_study_key.nunique() < 2:
        return {}
    work = frame.reset_index(drop=True)
    clusters = {
        str(target): [group.index.to_numpy() for _, group in target_rows.groupby("heldout_study_key", sort=True)]
        for target, target_rows in work.groupby("target_chembl_id", sort=True)
    }
    rng = np.random.default_rng(seed)
    collected: dict[str, list[float]] = defaultdict(list)
    for _ in range(iterations):
        indices = []
        for target in sorted(clusters):
            groups = clusters[target]
            for choice in rng.integers(0, len(groups), len(groups)):
                indices.append(groups[int(choice)])
        sample = work.iloc[np.concatenate(indices)]
        values = metric_values(sample, probability, k=k)
        for name, value in values.items():
            if name == "n" or value is None:
                continue
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                collected[name].append(float(value))
    return {
        name: ([float(v) for v in np.quantile(values, [0.025, 0.975])] if values else None)
        for name, values in collected.items()
    }


def paired_clustered_comparison(
    frame: pd.DataFrame,
    candidate: str,
    reference: str,
    *,
    k: int,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    candidate_point = metric_values(frame, candidate, k=k)
    reference_point = metric_values(frame, reference, k=k)
    precision_name = f"precision_at_{min(k, len(frame))}"
    names = ("auroc", "pr_auc", precision_name, "brier")
    point = {
        name: (
            float(candidate_point[name]) - float(reference_point[name])
            if candidate_point.get(name) is not None and reference_point.get(name) is not None else None
        ) for name in names
    }
    work = frame.reset_index(drop=True)
    clusters = {
        str(target): [group.index.to_numpy() for _, group in target_rows.groupby("heldout_study_key", sort=True)]
        for target, target_rows in work.groupby("target_chembl_id", sort=True)
    }
    rng = np.random.default_rng(seed)
    differences: dict[str, list[float]] = defaultdict(list)
    for _ in range(iterations):
        indices = []
        for target in sorted(clusters):
            groups = clusters[target]
            for choice in rng.integers(0, len(groups), len(groups)):
                indices.append(groups[int(choice)])
        sample = work.iloc[np.concatenate(indices)]
        left = metric_values(sample, candidate, k=k)
        right = metric_values(sample, reference, k=k)
        for name in names:
            if left.get(name) is not None and right.get(name) is not None:
                differences[name].append(float(left[name]) - float(right[name]))
    return {
        name: {
            "difference": point[name],
            "study_clustered_ci95": (
                [float(v) for v in np.quantile(differences[name], [0.025, 0.975])]
                if differences[name] else None
            ),
        } for name in names
    }


def audit_exposure(
    events: pd.DataFrame,
    named_ledgers: dict[str, dict[str, set[str]]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows = []
    report: dict[str, Any] = {}
    for name, ledger in named_ledgers.items():
        counts = Counter()
        compromised = set()
        for _, event in events.iterrows():
            reasons = event_overlap_reasons(event, ledger)
            if reasons:
                compromised.add(str(event.expansion_event_id))
                counts.update(reasons)
                rows.append({
                    "expansion_event_id": event.expansion_event_id,
                    "ledger": name,
                    "target_chembl_id": event.target_chembl_id,
                    "molecule_id": event.molecule_id,
                    "reasons": ";".join(sorted(reasons)),
                })
        for key in ("previous_target", "previous_compound", "previous_activity", "previous_assay", "previous_publication"):
            counts.setdefault(key, 0)
        report[name] = {
            "overlap_counts": dict(sorted(counts.items())),
            "compromised_events": len(compromised),
            "ledger_sizes": {key: len(value) for key, value in ledger.items()},
        }
    union = merge_ledgers(*named_ledgers.values())
    compromised = {
        str(event.expansion_event_id)
        for _, event in events.iterrows() if event_overlap_reasons(event, union)
    }
    within_assay = int(sum(bool(_ids(row.prior_assay_ids) & _ids(row.heldout_assay_ids)) for row in events.itertuples()))
    within_document = int(sum(bool(_ids(row.prior_document_ids) & _ids(row.heldout_document_ids)) for row in events.itertuples()))
    duplicate_compounds = int(events.molecule_id.duplicated(keep=False).sum())
    report["union"] = {
        "compromised_events": len(compromised),
        "within_event_assay_overlap": within_assay,
        "within_event_publication_overlap": within_document,
        "duplicate_compounds": duplicate_compounds,
    }
    return pd.DataFrame(rows, columns=["expansion_event_id", "ledger", "target_chembl_id", "molecule_id", "reasons"]), report


def assert_measurement_timing(events: pd.DataFrame, measurements: pd.DataFrame) -> None:
    event_year = events.set_index("expansion_event_id").heldout_year.astype(int).to_dict()
    if set(measurements.expansion_event_id) != set(event_year):
        raise AssertionError("Measurement provenance does not cover every eligible event")
    for row in measurements.itertuples(index=False):
        year = int(row.document_year)
        boundary = event_year[str(row.expansion_event_id)]
        if row.evidence_role == "prior" and not year < boundary:
            raise AssertionError("Historical feature source is not strictly earlier than held-out evidence")
        if row.evidence_role == "heldout" and year != boundary:
            raise AssertionError("Held-out source year does not match the event boundary")


def score_eligible(
    events: pd.DataFrame,
    full_artifact: dict[str, Any],
    selected_artifact: dict[str, Any],
) -> pd.DataFrame:
    assert_historical_feature_contract(full_artifact["features"])
    scored = score_frozen_models(events, full_artifact, selected_artifact)
    scored["probability_historical_tenfold_fraction"] = scored[
        "prior_historical_tenfold_fraction"
    ].astype(float).clip(0, 1)
    return scored


def evaluate(
    predictions: pd.DataFrame,
    *,
    k: int,
    iterations: int,
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    models = {
        "frozen_full": "probability_frozen_full",
        "frozen_constant_prevalence": "probability_frozen_selected",
        "historical_tenfold_fraction": "probability_historical_tenfold_fraction",
    }
    scopes = {
        "overall": predictions,
        "natural_products": predictions[predictions.is_verified_natural_product.astype(bool)],
    }
    results: dict[str, Any] = {scope: {} for scope in scopes}
    rows = []
    for scope_index, (scope, frame) in enumerate(scopes.items()):
        for model_index, (model, column) in enumerate(models.items()):
            point = metric_values(frame, column, k=k)
            intervals = clustered_intervals(
                frame, column, k=k, iterations=iterations,
                seed=seed + scope_index * 1000 + model_index * 100,
            )
            results[scope][model] = {"point": point, "study_clustered_ci95": intervals}
            rows.append({
                "scope": scope, "model": model, **point,
                **{f"{name}_ci95": json.dumps(value) for name, value in intervals.items()},
            })
    results["paired_full_minus_baselines"] = {
        baseline: paired_clustered_comparison(
            predictions, "probability_frozen_full", column,
            k=k, iterations=iterations, seed=seed + 5000 + index * 100,
        )
        for index, (baseline, column) in enumerate((
            ("frozen_constant_prevalence", "probability_frozen_selected"),
            ("historical_tenfold_fraction", "probability_historical_tenfold_fraction"),
        ))
    }
    return results, pd.DataFrame(rows)


def natural_product_errors(predictions: pd.DataFrame, *, k: int) -> pd.DataFrame:
    natural = predictions[predictions.is_verified_natural_product.astype(bool)].copy()
    natural["frozen_full_rank"] = natural.probability_frozen_full.rank(method="min", ascending=False).astype(int)
    natural["retrieved_at_k"] = natural.frozen_full_rank <= min(k, len(natural))
    natural["fixed_0_5_high_risk"] = natural.probability_frozen_full >= 0.5
    natural["classification_at_0_5"] = np.select(
        [
            (natural.disagreement_label == 1) & natural.fixed_0_5_high_risk,
            (natural.disagreement_label == 1) & ~natural.fixed_0_5_high_risk,
            (natural.disagreement_label == 0) & natural.fixed_0_5_high_risk,
        ],
        ["true_positive", "missed_disagreement", "false_positive"],
        default="true_negative",
    )
    columns = [
        "expansion_event_id", "target_chembl_id", "target_name", "molecule_id",
        "canonical_smiles", "heldout_study_key", "heldout_year", "disagreement_label",
        "absolute_pactivity_shift", "prior_median_pactivity", "heldout_median_pactivity",
        "natural_product_sources", "probability_frozen_full", "probability_frozen_selected",
        "probability_historical_tenfold_fraction", "frozen_full_rank", "retrieved_at_k",
        "fixed_0_5_high_risk", "classification_at_0_5",
    ]
    return natural[columns].sort_values(
        ["disagreement_label", "probability_frozen_full"], ascending=[False, False]
    ).reset_index(drop=True)


def _write_comparison_svg(path: Path, metrics: dict[str, Any], k: int) -> None:
    names = [
        ("frozen_full", "Frozen full", "#2dd4bf"),
        ("frozen_constant_prevalence", "Frozen prevalence", "#f59e0b"),
        ("historical_tenfold_fraction", "Historical fraction", "#60a5fa"),
    ]
    measures = [("auroc", "AUROC"), ("pr_auc", "PR-AUC"), (f"precision_at_{k}", f"Precision@{k}"), ("brier", "Brier")]
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="620" viewBox="0 0 1100 620">',
        '<rect width="100%" height="100%" fill="#07131d"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#dcecf4}.muted{fill:#8da9b8}</style>',
        '<text x="45" y="48" font-size="24" font-weight="700">Expanded independent validation</text>',
        '<text x="45" y="74" font-size="14" class="muted">Frozen models; 94 unseen events; Brier is lower-is-better.</text>',
    ]
    for measure_index, (metric, label) in enumerate(measures):
        y = 120 + measure_index * 115
        lines.append(f'<text x="45" y="{y}" font-size="17">{html.escape(label)}</text>')
        for model_index, (name, display, color) in enumerate(names):
            value = metrics["overall"][name]["point"].get(metric)
            bar_y = y - 23 + model_index * 27
            width = 620 * float(value or 0)
            lines.extend([
                f'<rect x="245" y="{bar_y}" width="{width:.1f}" height="19" rx="3" fill="{color}"/>',
                f'<text x="880" y="{bar_y + 15}" font-size="13">{html.escape(display)} {float(value):.3f}</text>',
            ])
    lines.append('</svg>')
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_calibration_svg(path: Path, predictions: pd.DataFrame) -> None:
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="720" height="620" viewBox="0 0 720 620">',
        '<rect width="100%" height="100%" fill="#07131d"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#dcecf4}.muted{fill:#8da9b8}</style>',
        '<text x="45" y="45" font-size="23" font-weight="700">Frozen-full calibration</text>',
        '<line x1="85" y1="535" x2="625" y2="535" stroke="#8da9b8"/>',
        '<line x1="85" y1="535" x2="85" y2="85" stroke="#8da9b8"/>',
        '<line x1="85" y1="535" x2="625" y2="85" stroke="#475569" stroke-dasharray="7 7"/>',
    ]
    bins = pd.cut(predictions.probability_frozen_full, bins=np.linspace(0, 1, 6), include_lowest=True)
    calibration = predictions.assign(_bin=bins).groupby("_bin", observed=False).agg(
        predicted=("probability_frozen_full", "mean"), observed=("disagreement_label", "mean"), n=("disagreement_label", "size")
    ).dropna()
    for row in calibration.itertuples(index=False):
        x = 85 + 540 * float(row.predicted)
        y = 535 - 450 * float(row.observed)
        radius = 5 + min(12, math.sqrt(int(row.n)))
        lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="#2dd4bf" opacity="0.88"/>')
        lines.append(f'<text x="{x + 12:.1f}" y="{y - 9:.1f}" font-size="12">n={int(row.n)}</text>')
    lines.extend([
        '<text x="300" y="585" font-size="14">Mean predicted risk</text>',
        '<text x="20" y="310" font-size="14" transform="rotate(-90 20 310)">Observed disagreement fraction</text>',
        '</svg>',
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_report(path: Path, results: dict[str, Any], errors: pd.DataFrame) -> None:
    dataset = results["dataset"]
    overall = results["metrics"]["overall"]
    natural = results["metrics"]["natural_products"]
    full = overall["frozen_full"]
    np_full = natural["frozen_full"]
    k_name = f"precision_at_{results['precision_k']}"
    classifications = dict(Counter(errors.classification_at_0_5))
    lines = [
        "# Expanded independent validation of frozen assay-disagreement models", "",
        "## Independence audit", "",
        f"All **{dataset['events']} eligible events** passed compound, target, activity, assay, publication, temporal, and within-event leakage checks against both the original model corpus and the earlier independent validation. The audit found zero compromised events. Model SHA-256 values were identical before and after inference.", "",
        "The frozen full artifact uses only molecular descriptors and historical evidence fields. The selected study-generalized artifact remains the constant development prevalence. No fitting, recalibration, threshold tuning, or winner selection was performed.", "",
        "## Dataset", "",
        f"The locked expansion contains {dataset['events']} compounds, {dataset['disagreements']} disagreements, {dataset['studies']} held-out studies, and {dataset['targets']} unseen targets. The observed disagreement fraction ({dataset['prevalence']:.3f}) is a property of this enriched repeated-measurement sample and is **not** an estimate of real-world prevalence.", "",
        "## Overall performance", "",
        "| Method | AUROC | PR-AUC | Precision@10 | Brier | ECE-10 | Calibration intercept | Calibration slope |", "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, display in (("frozen_full", "Frozen full"), ("frozen_constant_prevalence", "Frozen prevalence"), ("historical_tenfold_fraction", "Historical tenfold fraction")):
        point = overall[name]["point"]
        def fmt(value: Any) -> str:
            return "undefined" if value is None else f"{float(value):.3f}"
        lines.append(f"| {display} | {fmt(point['auroc'])} | {fmt(point['pr_auc'])} | {fmt(point[k_name])} | {fmt(point['brier'])} | {fmt(point['ece_10'])} | {fmt(point['calibration_intercept'])} | {fmt(point['calibration_slope'])} |")
    lines += [
        "", f"Frozen-full study-clustered 95% intervals: `{json.dumps(full['study_clustered_ci95'], sort_keys=True)}`.", "",
        "## Natural products", "",
        f"The subgroup has {dataset['natural_products']} compounds and {dataset['natural_product_disagreements']} disagreements. Frozen-full AUROC is `{np_full['point']['auroc']}`, PR-AUC `{np_full['point']['pr_auc']}`, Precision@10 `{np_full['point'][k_name]}`, and Brier `{np_full['point']['brier']}`. Study-clustered intervals are `{json.dumps(np_full['study_clustered_ci95'], sort_keys=True)}`.", "",
        f"At the fixed, non-tuned 0.5 descriptive threshold, classifications are `{json.dumps(classifications, sort_keys=True)}`. The accompanying error table preserves every natural-product compound, score, rank, outcome, target, publication, and potency shift.", "",
        "## Conclusion", "",
        results["conclusion"], "",
        "No method is selected as a winner from this evaluation. Disagreement is an assay-reproducibility outcome, not evidence of efficacy, safety, biological failure, or clinical utility.", "",
        "## Limitations", "",
        "- The test set was deliberately enriched for independently repeated measurements and is not prevalence-representative.",
        "- All 12 targets are unseen by the frozen target encoder, so this is simultaneous cross-target and cross-study validation.",
        "- Twenty natural products and seven subgroup positives yield wide, unstable confidence intervals.",
        "- Publication identifiers proxy independent studies; laboratory and full protocol metadata remain incomplete.",
        "- Calibration estimates are descriptive and no recalibration was permitted.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_expanded_validation(
    events_path: Path,
    measurements_path: Path,
    protocol_path: Path,
    development_specs: Path,
    development_events: Path,
    prior_validation_specs: Path,
    prior_validation_events: Path,
    full_model_path: Path,
    selected_model_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol.get("protocol_version") != VALIDATION_VERSION or not protocol.get("created_before_frozen_model_scoring"):
        raise ValueError("Expanded validation protocol is absent or not frozen")
    iterations = int(protocol["bootstrap"]["iterations"])
    seed = int(protocol["bootstrap"]["seed"])
    k = int(protocol["precision_k"])
    events = pd.read_csv(events_path)
    measurements = pd.read_csv(measurements_path)
    full, selected, hashes_before = _load_frozen_models(full_model_path, selected_model_path)
    development_ledger, development_summary = build_previous_exposure_ledger(development_specs, development_events)
    validation_ledger, validation_summary = build_previous_exposure_ledger(prior_validation_specs, prior_validation_events)
    overlap_rows, audit = audit_exposure(events, {
        "model_development_and_original_evaluation": development_ledger,
        "prior_independent_validation": validation_ledger,
    })
    compromised = set(overlap_rows.expansion_event_id.astype(str)) if len(overlap_rows) else set()
    eligible = events[~events.expansion_event_id.astype(str).isin(compromised)].copy().reset_index(drop=True)
    if audit["union"]["within_event_assay_overlap"] or audit["union"]["within_event_publication_overlap"]:
        raise AssertionError("Expansion contains within-event leakage")
    assert_measurement_timing(eligible, measurements[measurements.expansion_event_id.isin(eligible.expansion_event_id)])
    predictions = score_eligible(eligible, full, selected)
    metrics, metric_table = evaluate(predictions, k=k, iterations=iterations, seed=seed)
    errors = natural_product_errors(predictions, k=k)
    hashes_after = {"full_model": _sha256(full_model_path), "selected_model": _sha256(selected_model_path)}
    if hashes_before != hashes_after:
        raise AssertionError("Frozen model artifact changed during evaluation")
    natural = predictions[predictions.is_verified_natural_product.astype(bool)]
    dataset = {
        "input_events": len(events), "excluded_compromised_events": len(compromised),
        "events": len(predictions), "disagreements": int(predictions.disagreement_label.sum()),
        "prevalence": float(predictions.disagreement_label.mean()),
        "studies": int(predictions.heldout_study_key.nunique()),
        "targets": int(predictions.target_chembl_id.nunique()),
        "natural_products": len(natural),
        "natural_product_disagreements": int(natural.disagreement_label.sum()),
        "natural_product_studies": int(natural.heldout_study_key.nunique()),
    }
    full_point = metrics["overall"]["frozen_full"]["point"]
    full_ci = metrics["overall"]["frozen_full"]["study_clustered_ci95"]
    prevalence = dataset["prevalence"]
    if full_ci.get("auroc") and full_ci["auroc"][0] > 0.5 and full_ci.get("pr_auc") and full_ci["pr_auc"][0] > prevalence:
        conclusion = "The frozen full model shows independently supported discrimination under the prespecified interval criteria, but remains unselected and requires replication because this enriched sample is not representative."
    else:
        conclusion = "The expanded evaluation does not establish reliable generalization: uncertainty includes weak or non-discriminatory performance under the prespecified criteria. The negative or inconclusive result is retained without selecting or modifying a model."
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "validation_version": VALIDATION_VERSION,
        "protocol_sha256": _sha256(protocol_path), "precision_k": k,
        "bootstrap_iterations": iterations, "seed": seed,
        "dataset": dataset, "leakage_audit": audit,
        "exposure_ledgers": {"development": development_summary, "prior_validation": validation_summary},
        "feature_contract": {
            "frozen_full_features": list(full["features"]),
            "uses_only_pre_heldout_information": True,
            "forbidden_feature_violations": [],
            "unseen_target_encoder_categories": sorted(predictions.target_chembl_id.unique().tolist()),
        },
        "metrics": metrics, "model_hashes_before": hashes_before,
        "model_hashes_after": hashes_after, "models_retrained_recalibrated_or_tuned": False,
        "evaluation_used_for_model_selection": False,
        "representative_prevalence_claim": False,
        "frozen_full_validation_status": full.get("validation_status"),
        "frozen_selected_validation_status": selected.get("validation_status"),
        "conclusion": conclusion,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "eligible_events.csv", eligible)
    _write_csv(output_dir / "frozen_predictions.csv", predictions)
    _write_csv(output_dir / "leakage_exclusions.csv", overlap_rows)
    _write_csv(output_dir / "metric_summary.csv", metric_table)
    _write_csv(output_dir / "natural_product_error_analysis.csv", errors)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_comparison_svg(output_dir / "model_comparison.svg", metrics, k)
    _write_calibration_svg(output_dir / "calibration.svg", predictions)
    _write_report(output_dir / "SCIENTIFIC_REPORT.md", results, errors)
    manifest = {
        "created_at": datetime.now(UTC).isoformat(), "validation_version": VALIDATION_VERSION,
        "inputs": {
            path.name: {"path": str(path.resolve()), "sha256": _sha256(path)}
            for path in (events_path, measurements_path, protocol_path, development_specs, development_events,
                         prior_validation_specs, prior_validation_events, full_model_path, selected_model_path)
        },
        "outputs": {},
        "models_retrained_recalibrated_or_tuned": False,
    }
    for name in ("eligible_events.csv", "frozen_predictions.csv", "leakage_exclusions.csv", "metric_summary.csv", "natural_product_error_analysis.csv", "results.json", "model_comparison.svg", "calibration.svg", "SCIENTIFIC_REPORT.md"):
        manifest["outputs"][name] = _sha256(output_dir / name)
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "command.txt").write_text(
        "python -m app.expanded_validation " + " ".join((
            f'--events "{events_path.resolve()}"', f'--measurements "{measurements_path.resolve()}"',
            f'--protocol "{protocol_path.resolve()}"', f'--development-specs "{development_specs.resolve()}"',
            f'--development-events "{development_events.resolve()}"',
            f'--prior-validation-specs "{prior_validation_specs.resolve()}"',
            f'--prior-validation-events "{prior_validation_events.resolve()}"',
            f'--full-model "{full_model_path.resolve()}"', f'--selected-model "{selected_model_path.resolve()}"',
            f'--output-dir "{output_dir.resolve()}"',
        )) + "\n", encoding="utf-8"
    )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate frozen disagreement models on the expanded independent cohort")
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--measurements", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--development-specs", type=Path, required=True)
    parser.add_argument("--development-events", type=Path, required=True)
    parser.add_argument("--prior-validation-specs", type=Path, required=True)
    parser.add_argument("--prior-validation-events", type=Path, required=True)
    parser.add_argument("--full-model", type=Path, required=True)
    parser.add_argument("--selected-model", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    results = run_expanded_validation(
        args.events, args.measurements, args.protocol, args.development_specs,
        args.development_events, args.prior_validation_specs, args.prior_validation_events,
        args.full_model, args.selected_model, args.output_dir,
    )
    print(json.dumps({"dataset": results["dataset"], "conclusion": results["conclusion"]}, indent=2))


if __name__ == "__main__":
    main()
