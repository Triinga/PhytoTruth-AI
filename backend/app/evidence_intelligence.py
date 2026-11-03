from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


SCHEMES = ("study_grouped", "forward_temporal")
METHODS = (
    "median_potency",
    "student_t_regression",
    "xgboost",
    "phytotrust_m0",
    "m1_no_protocol",
    "phytotrust_m1",
)
THRESHOLD = 1.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_robust_method(metrics: dict[str, Any]) -> dict[str, Any]:
    """Minimize the worst NLL rank across two prespecified validation schemes."""
    scheme_ranks: dict[str, dict[str, int]] = {}
    for scheme in SCHEMES:
        available = metrics[scheme]["all"]
        ordered = sorted(METHODS, key=lambda method: (available[method]["nll"], method))
        scheme_ranks[scheme] = {method: rank for rank, method in enumerate(ordered, 1)}
    candidates = []
    for method in METHODS:
        ranks = {scheme: scheme_ranks[scheme][method] for scheme in SCHEMES}
        calibration = float(np.mean([
            metrics[scheme]["all"][method]["mean_absolute_calibration_error"]
            for scheme in SCHEMES
        ]))
        candidates.append({
            "method": method,
            "ranks": ranks,
            "worst_nll_rank": max(ranks.values()),
            "mean_nll_rank": float(np.mean(list(ranks.values()))),
            "mean_calibration_error": calibration,
        })
    selected = min(
        candidates,
        key=lambda row: (
            row["worst_nll_rank"], row["mean_nll_rank"], row["mean_calibration_error"], row["method"]
        ),
    )
    return {
        "rule": "Minimize worst NLL rank across study-grouped and forward-temporal validation; break ties by mean NLL rank, then calibration error.",
        "selected_method": selected["method"],
        "selected": selected,
        "candidates": candidates,
    }


def disagreement_probability(row: pd.Series, *, threshold: float = THRESHOLD, df: float = 4.0) -> float:
    """Probability that a future value differs from the historical median by threshold."""
    center = float(row.prediction)
    historical = float(row.prior_median_pactivity)
    half_width = max(float(row.upper_90) - float(row.lower_90), 1e-9) / 2
    scale = max(half_width / student_t.ppf(0.95, df), 1e-6)
    lower = historical - threshold
    upper = historical + threshold
    return float(student_t.cdf(lower, df=df, loc=center, scale=scale) + student_t.sf(upper, df=df, loc=center, scale=scale))


def expected_information_gain(
    interval_width_90: float,
    prior_study_count: int,
    *,
    df: float = 4.0,
    observation_noise_sd: float = 0.50,
) -> float:
    """Gaussian-approximate entropy reduction from one independent observation."""
    scale = max(interval_width_90 / (2 * student_t.ppf(0.95, df)), 1e-9)
    predictive_variance = scale * scale * df / (df - 2)
    latent_variance = predictive_variance / max(prior_study_count, 1)
    noise_variance = observation_noise_sd * observation_noise_sd
    return float(0.5 * math.log1p(latent_variance / noise_variance))


def _cluster_bootstrap_ci(
    frame: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], float],
    *,
    iterations: int,
    seed: int,
) -> list[float]:
    studies = sorted(frame.heldout_study_key.astype(str).unique())
    groups = {study: frame[frame.heldout_study_key.astype(str) == study] for study in studies}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(iterations):
        sampled = rng.choice(studies, size=len(studies), replace=True)
        boot = pd.concat([groups[study] for study in sampled], ignore_index=True)
        try:
            value = float(statistic(boot))
            if math.isfinite(value):
                values.append(value)
        except ValueError:
            continue
    if not values:
        return [math.nan, math.nan]
    return [float(value) for value in np.quantile(values, [0.025, 0.975])]


def _binary_metrics(frame: pd.DataFrame, *, iterations: int, seed: int) -> dict[str, Any]:
    labels = frame.disagreement_label.astype(int)
    probabilities = frame.disagreement_probability.astype(float)
    def _auroc(data: pd.DataFrame) -> float:
        if data.disagreement_label.nunique() < 2:
            raise ValueError("AUROC requires both classes")
        return float(roc_auc_score(data.disagreement_label, data.disagreement_probability))

    def _average_precision(data: pd.DataFrame) -> float:
        if data.disagreement_label.nunique() < 2:
            raise ValueError("Average precision requires both classes")
        return float(average_precision_score(data.disagreement_label, data.disagreement_probability))

    metric_functions: dict[str, Callable[[pd.DataFrame], float]] = {
        "brier": lambda data: brier_score_loss(data.disagreement_label, data.disagreement_probability),
        "auroc": _auroc,
        "average_precision": _average_precision,
    }
    output: dict[str, Any] = {
        "n": int(len(frame)),
        "studies": int(frame.heldout_study_key.nunique()),
        "disagreements": int(labels.sum()),
        "prevalence": float(labels.mean()),
    }
    for index, (name, function) in enumerate(metric_functions.items()):
        try:
            point = float(function(frame))
        except ValueError:
            point = math.nan
        output[name] = point
        output[f"{name}_study_cluster_ci95"] = _cluster_bootstrap_ci(
            frame, function, iterations=iterations, seed=seed + index
        )
    return output


def _selective_curve(frame: pd.DataFrame) -> list[dict[str, Any]]:
    ordered = frame.sort_values("interval_width_90", kind="stable")
    rows = []
    for retained_fraction in (1.0, 0.8, 0.6, 0.4):
        retained = ordered.iloc[: max(1, math.ceil(len(ordered) * retained_fraction))]
        rows.append({
            "retained_fraction": retained_fraction,
            "retained_n": int(len(retained)),
            "abstained_n": int(len(frame) - len(retained)),
            "mae": float(np.mean(np.abs(retained.prediction - retained.heldout_median_pactivity))),
            "coverage_90": float(np.mean(
                (retained.heldout_median_pactivity >= retained.lower_90)
                & (retained.heldout_median_pactivity <= retained.upper_90)
            )),
            "mean_interval_width_90": float(retained.interval_width_90.mean()),
        })
    return rows


def _recommendation(row: pd.Series) -> dict[str, Any]:
    probability = float(row.disagreement_probability)
    missingness = float(row.protocol_missing_fraction)
    if missingness >= 0.9:
        experiment = "Protocol-complete independent IC50 replication"
        rationale = "Historical protocol metadata are mostly missing; pre-register construct, substrate, detection method, duration, temperature, and pH."
    elif probability >= 0.35:
        experiment = "Orthogonal-method independent IC50 replication"
        rationale = "The predictive distribution assigns material probability to a ten-fold or larger future-study shift."
    else:
        experiment = "Independent confirmatory IC50 replication"
        rationale = "A new publication-level replicate would reduce uncertainty without reusing the historical assay batch."
    return {
        "target_chembl_id": str(row.target_chembl_id),
        "molecule_id": str(row.molecule_id),
        "prior_median_pactivity": float(row.prior_median_pactivity),
        "predicted_future_pactivity": float(row.prediction),
        "prediction_interval_90": [float(row.lower_90), float(row.upper_90)],
        "disagreement_probability": probability,
        "prior_studies": int(row.prior_study_count),
        "protocol_missing_fraction": missingness,
        "expected_information_gain_nats": float(row.expected_information_gain_nats),
        "recommended_experiment": experiment,
        "rationale": rationale,
        "decision_status": "research_planning_only_not_validated_clinical_or_biological_truth",
    }


def compile_intelligence(
    source_results: Path,
    predictions_path: Path,
    output_dir: Path,
    *,
    bootstrap_iterations: int = 5_000,
    seed: int = 2026,
) -> dict[str, Any]:
    source = json.loads(source_results.read_text(encoding="utf-8"))
    predictions = pd.read_csv(predictions_path)
    selection = select_robust_method(source["metrics"])
    selected_method = selection["selected_method"]
    evaluated = predictions[
        (predictions.validation_scheme == "forward_temporal")
        & (predictions.method == selected_method)
    ].copy()
    if evaluated.empty:
        raise ValueError(f"No forward-temporal predictions for {selected_method}")
    evaluated["interval_width_90"] = evaluated.upper_90 - evaluated.lower_90
    evaluated["disagreement_label"] = (
        np.abs(evaluated.heldout_median_pactivity - evaluated.prior_median_pactivity) >= THRESHOLD
    ).astype(int)
    evaluated["disagreement_probability"] = evaluated.apply(disagreement_probability, axis=1)
    evaluated["expected_information_gain_nats"] = evaluated.apply(
        lambda row: expected_information_gain(row.interval_width_90, int(row.prior_study_count)), axis=1
    )
    binary = _binary_metrics(evaluated, iterations=bootstrap_iterations, seed=seed)
    selective = _selective_curve(evaluated)

    # The prioritization key is deliberately outcome-blind: no held-out value or label is read here.
    prioritized = evaluated.sort_values(
        ["expected_information_gain_nats", "disagreement_probability", "molecule_id"],
        ascending=[False, False, True],
        kind="stable",
    ).drop_duplicates("molecule_id")
    recommendations = [_recommendation(row) for _, row in prioritized.head(8).iterrows()]
    result = {
        "experiment": "PhytoTruth Study-Aware Evidence Intelligence",
        "generated_at": datetime.now(UTC).isoformat(),
        "seed": seed,
        "source_model_version": source["model_version"],
        "source_status": source["scientific_status"],
        "selection": selection,
        "forward_temporal_distribution_metrics": source["metrics"]["forward_temporal"]["all"][selected_method],
        "study_grouped_distribution_metrics": source["metrics"]["study_grouped"]["all"][selected_method],
        "disagreement_evaluation": binary,
        "selective_prediction": selective,
        "recommendations": recommendations,
        "integrity": {
            "heldout_outcome_used_for_recommendation_priority": False,
            "natural_product_status_used_as_predictor": False,
            "heldout_protocol_used_as_predictor": False,
            "confirmatory_status": "retrospective_development_only_secondary_analysis",
            "deployment_status": "research_planning_guardrail_not_validated_for_autonomous_decisions",
            "warning": "Expected information gain is a model-based Gaussian approximation. It prioritizes information collection; it is not a probability of efficacy, safety, or experimental success.",
        },
        "scientific_conclusion": (
            "The robust probabilistic guardrail supplies calibrated-looking intervals and explicit abstention, but this secondary retrospective analysis does not establish independent generalization. Recommendations are auditable research-planning hypotheses, not validated biological decisions."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    evaluated.to_csv(output_dir / "evaluated_predictions.csv", index=False)
    (output_dir / "recommendations.json").write_text(json.dumps(recommendations, indent=2), encoding="utf-8")
    manifest = {
        "source_results": str(source_results.resolve()),
        "source_results_sha256": _sha256(source_results),
        "predictions": str(predictions_path.resolve()),
        "predictions_sha256": _sha256(predictions_path),
        "results_sha256": _sha256(output_dir / "results.json"),
        "outcome_blind_recommendation_fields": [
            "prediction", "lower_90", "upper_90", "prior_median_pactivity", "prior_study_count", "protocol_missing_fraction"
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    (output_dir / "command.txt").write_text(
        "python -m app.evidence_intelligence "
        f"--source-results {source_results} --predictions {predictions_path} "
        f"--output-dir {output_dir} --bootstrap-iterations {bootstrap_iterations} --seed {seed}\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile study-aware evidence intelligence artifacts")
    parser.add_argument("--source-results", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    result = compile_intelligence(
        args.source_results,
        args.predictions,
        args.output_dir,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    print(json.dumps({
        "selected_method": result["selection"]["selected_method"],
        "disagreement_evaluation": result["disagreement_evaluation"],
        "status": result["integrity"]["deployment_status"],
    }, indent=2))


if __name__ == "__main__":
    main()
