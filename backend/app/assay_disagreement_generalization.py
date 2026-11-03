from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import joblib
import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .assay_disagreement import (
    ACTIVE_THRESHOLD,
    ALL_NUMERIC_FEATURES,
    DEVELOPMENT_END_YEAR,
    EVALUATION_START_YEAR,
    MODEL_FEATURES,
    OUTCOME_THRESHOLD_LOG10,
    _metric_values,
    _sha256,
    _write_csv,
    build_dataset,
    clustered_metric_intervals,
    clustered_model_comparison,
)


ENHANCED_MOLECULAR = [
    "molecular_weight", "logp", "tpsa", "fraction_csp3", "ring_count",
    "aromatic_atom_fraction", "heteroatom_fraction", "bertz_complexity",
]
COMPACT_HISTORY = [
    "prior_median_pactivity", "prior_distance_to_active_threshold", "prior_documents",
    "prior_assays", "prior_year_span", "prior_mad", "prior_document_median_range",
    "prior_assay_median_range", "prior_historical_tenfold_fraction",
    "prior_recent_shift", "prior_temporal_slope",
]
COMPACT_METADATA = [
    "prior_records_per_document", "prior_assays_per_document", "prior_journals",
    "prior_descriptions", "prior_missing_description_fraction", "prior_missing_journal_fraction",
    "prior_cell_keyword_fraction", "prior_recombinant_keyword_fraction",
    "prior_substrate_keyword_fraction",
]
COMPACT_FEATURES = ENHANCED_MOLECULAR + COMPACT_HISTORY + COMPACT_METADATA

CANDIDATES: dict[str, dict[str, Any]] = {
    "current_full_logistic": {"kind": "logistic", "features": list(ALL_NUMERIC_FEATURES), "calibrate": False},
    "current_full_calibrated": {"kind": "logistic", "features": list(ALL_NUMERIC_FEATURES), "calibrate": True},
    "compact_history_logistic": {"kind": "logistic", "features": list(COMPACT_HISTORY + COMPACT_METADATA), "calibrate": True},
    "compact_logistic": {"kind": "logistic", "features": list(COMPACT_FEATURES), "calibrate": True},
    "compact_hist_gradient_boosting": {"kind": "hist_gradient", "features": list(COMPACT_FEATURES), "calibrate": True},
}
SELECTABLE_MODELS = (
    "constant_prevalence", "current_full_calibrated", "compact_history_logistic",
    "compact_logistic", "compact_hist_gradient_boosting",
)


def _model_pipeline(kind: str, features: Sequence[str], seed: int) -> Pipeline:
    numeric = Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())])
    if kind == "logistic":
        target = OneHotEncoder(handle_unknown="ignore")
        estimator = LogisticRegression(C=1.0, solver="liblinear", max_iter=2000, random_state=seed)
    elif kind == "hist_gradient":
        target = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
        estimator = HistGradientBoostingClassifier(
            learning_rate=0.05, max_iter=150, max_leaf_nodes=7, min_samples_leaf=15,
            l2_regularization=1.0, random_state=seed,
        )
    else:
        raise ValueError(f"Unknown model kind: {kind}")
    return Pipeline([
        ("preprocess", ColumnTransformer([
            ("numeric", numeric, list(features)),
            ("target", target, ["target_chembl_id"]),
        ])),
        ("classifier", estimator),
    ])


def _columns(features: Sequence[str]) -> list[str]:
    return list(features) + ["target_chembl_id"]


def leave_one_study_out_indices(studies: Sequence[str]) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = np.asarray(studies, dtype=str)
    return [
        (np.flatnonzero(groups != study), np.flatnonzero(groups == study))
        for study in sorted(set(groups))
    ]


def leave_one_study_out_predictions(
    development: pd.DataFrame,
    *,
    kind: str,
    features: Sequence[str],
    seed: int,
) -> np.ndarray:
    """Every prediction comes from a model that never saw that publication group."""
    probabilities = np.full(len(development), np.nan, dtype=float)
    groups = development["heldout_study_key"].astype(str).to_numpy()
    for fold, (train_indices, test_indices) in enumerate(leave_one_study_out_indices(groups)):
        training = development.iloc[train_indices]
        if training.disagreement_label.nunique() < 2:
            probabilities[test_indices] = float(training.disagreement_label.mean())
            continue
        model = _model_pipeline(kind, features, seed + fold)
        model.fit(training[_columns(features)], training.disagreement_label)
        probabilities[test_indices] = model.predict_proba(
            development.iloc[test_indices][_columns(features)]
        )[:, 1]
    if not np.isfinite(probabilities).all():
        raise AssertionError("LOSO prediction coverage is incomplete")
    return probabilities


def cross_fitted_calibration(
    raw_probabilities: Sequence[float], labels: Sequence[int], studies: Sequence[str], *, seed: int
) -> tuple[np.ndarray, LogisticRegression]:
    """Calibrate LOSO predictions without fitting a calibrator on its own study."""
    raw = np.clip(np.asarray(raw_probabilities, dtype=float), 1e-6, 1 - 1e-6)
    y = np.asarray(labels, dtype=int)
    groups = np.asarray(studies, dtype=str)
    splitter = GroupKFold(n_splits=min(5, len(set(groups))))
    calibrated = np.full(len(y), np.nan, dtype=float)
    for train_indices, test_indices in splitter.split(raw, y, groups):
        calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000, random_state=seed)
        calibrator.fit(logit(raw[train_indices]).reshape(-1, 1), y[train_indices])
        calibrated[test_indices] = calibrator.predict_proba(logit(raw[test_indices]).reshape(-1, 1))[:, 1]
    final_calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000, random_state=seed)
    final_calibrator.fit(logit(raw).reshape(-1, 1), y)
    if not np.isfinite(calibrated).all():
        raise AssertionError("Cross-fitted calibration coverage is incomplete")
    return calibrated, final_calibrator


def study_balanced_metrics(
    labels: Sequence[int], probabilities: Sequence[float], studies: Sequence[str]
) -> dict[str, float | None]:
    frame = pd.DataFrame({"label": labels, "probability": probabilities, "study": studies})
    sizes = frame.groupby("study").size()
    weights = frame.study.map(lambda study: 1.0 / sizes.loc[study]).to_numpy(float).copy()
    weights /= weights.sum()
    y = frame.label.to_numpy(int)
    p = np.clip(frame.probability.to_numpy(float), 1e-6, 1 - 1e-6)
    return {
        "n": len(frame),
        "studies": len(sizes),
        "auroc": float(roc_auc_score(y, p, sample_weight=weights)) if len(set(y)) == 2 else None,
        "pr_auc": float(average_precision_score(y, p, sample_weight=weights)) if len(set(y)) == 2 else None,
        "brier": float(np.sum(weights * np.square(y - p))),
        "log_loss": float(log_loss(y, p, sample_weight=weights, labels=[0, 1])),
    }


def development_model_selection(
    development: pd.DataFrame, *, seed: int
) -> tuple[dict[str, Any], pd.DataFrame, dict[str, LogisticRegression | None]]:
    oof = development[[
        "molecule_id", "target_chembl_id", "heldout_study_key", "heldout_year",
        "disagreement_label", "is_verified_natural_product",
    ]].copy()
    calibrators: dict[str, LogisticRegression | None] = {}
    prevalence = float(development.disagreement_label.mean())
    oof["probability_constant_prevalence"] = prevalence
    calibrators["constant_prevalence"] = None
    metrics: dict[str, Any] = {
        "constant_prevalence": {
            "pooled": _metric_values(development.disagreement_label, oof.probability_constant_prevalence),
            "study_balanced": study_balanced_metrics(
                development.disagreement_label, oof.probability_constant_prevalence,
                development.heldout_study_key,
            ),
        }
    }
    for index, (name, config) in enumerate(CANDIDATES.items()):
        raw = leave_one_study_out_predictions(
            development, kind=config["kind"], features=config["features"], seed=seed + index * 1000
        )
        oof[f"raw_probability_{name}"] = raw
        if config["calibrate"]:
            probabilities, calibrator = cross_fitted_calibration(
                raw, development.disagreement_label, development.heldout_study_key,
                seed=seed + index,
            )
        else:
            probabilities, calibrator = raw, None
        oof[f"probability_{name}"] = probabilities
        calibrators[name] = calibrator
        metrics[name] = {
            "pooled": _metric_values(development.disagreement_label, probabilities),
            "study_balanced": study_balanced_metrics(
                development.disagreement_label, probabilities, development.heldout_study_key,
            ),
        }
    selected = min(
        SELECTABLE_MODELS,
        key=lambda name: (float(metrics[name]["study_balanced"]["brier"]), name),
    )
    return {
        "selected_model": selected,
        "criterion": "minimum equal-study-weighted Brier score from leave-one-study-out development predictions",
        "candidate_metrics": metrics,
        "evaluation_used_for_selection": False,
    }, oof, calibrators


def fit_and_predict_evaluation(
    development: pd.DataFrame,
    evaluation: pd.DataFrame,
    selection: dict[str, Any],
    calibrators: dict[str, LogisticRegression | None],
    *,
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    predictions = evaluation.copy()
    fitted: dict[str, Any] = {"constant_prevalence": float(development.disagreement_label.mean())}
    predictions["probability_constant_prevalence"] = fitted["constant_prevalence"]
    for index, (name, config) in enumerate(CANDIDATES.items()):
        model = _model_pipeline(config["kind"], config["features"], seed + index * 1000)
        model.fit(development[_columns(config["features"])], development.disagreement_label)
        raw = model.predict_proba(evaluation[_columns(config["features"])])[:, 1]
        predictions[f"raw_probability_{name}"] = raw
        calibrator = calibrators[name]
        probabilities = (
            calibrator.predict_proba(logit(np.clip(raw, 1e-6, 1 - 1e-6)).reshape(-1, 1))[:, 1]
            if calibrator is not None else raw
        )
        predictions[f"probability_{name}"] = probabilities
        fitted[name] = {"model": model, "calibrator": calibrator, "config": config}
    predictions["selected_probability"] = predictions[f"probability_{selection['selected_model']}"]
    return fitted, predictions


def _study_table(frame: pd.DataFrame, probability_column: str) -> list[dict[str, Any]]:
    overall = _metric_values(frame.disagreement_label, frame[probability_column])
    rows: list[dict[str, Any]] = []
    for study, group in frame.groupby("heldout_study_key", sort=True):
        without = frame[frame.heldout_study_key != study]
        without_metrics = _metric_values(without.disagreement_label, without[probability_column])
        rows.append({
            "study": study,
            "target_chembl_id": ";".join(sorted(group.target_chembl_id.unique())),
            "year": int(group.heldout_year.min()),
            "molecules": len(group),
            "disagreements": int(group.disagreement_label.sum()),
            "natural_products": int(group.is_verified_natural_product.sum()),
            "observed_rate": float(group.disagreement_label.mean()),
            "mean_predicted_risk": float(group[probability_column].mean()),
            "brier": float(brier_score_loss(group.disagreement_label, group[probability_column])),
            "auroc_without_study": without_metrics["auroc"],
            "auroc_influence": (
                float(overall["auroc"]) - float(without_metrics["auroc"])
                if overall["auroc"] is not None and without_metrics["auroc"] is not None else None
            ),
            "pr_auc_without_study": without_metrics["pr_auc"],
            "pr_auc_influence": (
                float(overall["pr_auc"]) - float(without_metrics["pr_auc"])
                if overall["pr_auc"] is not None and without_metrics["pr_auc"] is not None else None
            ),
        })
    return rows


def _binary_icc(frame: pd.DataFrame) -> float | None:
    groups = [group.disagreement_label.to_numpy(float) for _, group in frame.groupby("heldout_study_key")]
    if len(groups) < 2 or len(frame) <= len(groups):
        return None
    grand = float(frame.disagreement_label.mean())
    ss_between = sum(len(values) * (float(values.mean()) - grand) ** 2 for values in groups)
    ss_within = sum(float(np.square(values - values.mean()).sum()) for values in groups)
    ms_between = ss_between / (len(groups) - 1)
    ms_within = ss_within / (len(frame) - len(groups))
    k_effective = (len(frame) - sum(len(values) ** 2 for values in groups) / len(frame)) / (len(groups) - 1)
    denominator = ms_between + (k_effective - 1) * ms_within
    return float((ms_between - ms_within) / denominator) if denominator else None


def _applicability_domain(development: pd.DataFrame, evaluation: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, float]:
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    train = scaler.fit_transform(imputer.fit_transform(development[ENHANCED_MOLECULAR]))
    test = scaler.transform(imputer.transform(evaluation[ENHANCED_MOLECULAR]))
    neighbors = NearestNeighbors(n_neighbors=2).fit(train)
    train_distance = neighbors.kneighbors(train)[0][:, 1]
    threshold = float(np.quantile(train_distance, 0.95))
    test_distance = NearestNeighbors(n_neighbors=1).fit(train).kneighbors(test)[0][:, 0]
    return test_distance, test_distance > threshold, threshold


def _natural_product_analysis(
    development: pd.DataFrame, predictions: pd.DataFrame
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    distances, outside, threshold = _applicability_domain(development, predictions)
    predictions = predictions.copy()
    predictions["descriptor_nn_distance"] = distances
    predictions["outside_descriptor_domain"] = outside
    probabilities = predictions.selected_probability.to_numpy(float)
    labels = predictions.disagreement_label.to_numpy(int)
    predicted = probabilities >= 0.5
    predictions["prediction_error_type"] = np.select(
        [predicted & (labels == 1), predicted & (labels == 0), ~predicted & (labels == 1)],
        ["true_positive", "false_positive", "false_negative"],
        default="true_negative",
    )
    predictions["squared_error"] = np.square(labels - probabilities)
    current_probabilities = predictions.probability_current_full_logistic.to_numpy(float)
    current_predicted = current_probabilities >= 0.5
    predictions["current_full_error_type"] = np.select(
        [current_predicted & (labels == 1), current_predicted & (labels == 0), ~current_predicted & (labels == 1)],
        ["true_positive", "false_positive", "false_negative"],
        default="true_negative",
    )
    predictions["current_full_squared_error"] = np.square(labels - current_probabilities)
    natural = predictions[predictions.is_verified_natural_product].copy()
    shift_rows = []
    for feature in COMPACT_FEATURES:
        left = predictions.loc[predictions.is_verified_natural_product, feature].astype(float)
        right = predictions.loc[~predictions.is_verified_natural_product, feature].astype(float)
        pooled = float(predictions[feature].astype(float).std(ddof=0))
        shift_rows.append({
            "feature": feature,
            "natural_product_mean": float(left.mean()),
            "other_mean": float(right.mean()),
            "standardized_mean_difference": float((left.mean() - right.mean()) / pooled) if pooled else 0.0,
        })
    summary = {
        "n": len(natural),
        "positives": int(natural.disagreement_label.sum()),
        "error_types": {str(key): int(value) for key, value in natural.prediction_error_type.value_counts().items()},
        "current_full_error_types": {
            str(key): int(value) for key, value in natural.current_full_error_type.value_counts().items()
        },
        "mean_squared_error": float(natural.squared_error.mean()),
        "current_full_mean_squared_error": float(natural.current_full_squared_error.mean()),
        "outside_descriptor_domain": int(natural.outside_descriptor_domain.sum()),
        "outside_descriptor_domain_fraction": float(natural.outside_descriptor_domain.mean()),
        "descriptor_domain_threshold": threshold,
        "by_target": {
            str(target): {
                "n": len(group), "positives": int(group.disagreement_label.sum()),
                "brier": float(group.squared_error.mean()),
            }
            for target, group in natural.groupby("target_chembl_id", sort=True)
        },
    }
    return natural, shift_rows, summary


def _bootstrap_selected_intervals(
    development: pd.DataFrame,
    evaluation: pd.DataFrame,
    model_name: str,
    calibrator: LogisticRegression | None,
    *,
    iterations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    if model_name == "constant_prevalence":
        estimates = np.full((iterations, len(evaluation)), float(development.disagreement_label.mean()))
        return estimates[0], estimates[0]
    config = CANDIDATES[model_name]
    rng = np.random.default_rng(seed)
    groups = [group.index.to_numpy() for _, group in development.groupby("heldout_study_key", sort=True)]
    estimates = []
    for iteration in range(iterations):
        choices = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[int(choice)] for choice in choices])
        sample = development.loc[indices]
        if sample.disagreement_label.nunique() < 2:
            continue
        model = _model_pipeline(config["kind"], config["features"], seed + iteration)
        model.fit(sample[_columns(config["features"])], sample.disagreement_label)
        raw = model.predict_proba(evaluation[_columns(config["features"])])[:, 1]
        estimate = (
            calibrator.predict_proba(logit(np.clip(raw, 1e-6, 1 - 1e-6)).reshape(-1, 1))[:, 1]
            if calibrator is not None else raw
        )
        estimates.append(estimate)
    matrix = np.asarray(estimates)
    return np.quantile(matrix, 0.025, axis=0), np.quantile(matrix, 0.975, axis=0)


def _write_plot(path: Path, results: dict[str, Any]) -> None:
    names = ["current_full_logistic", results["development_selection"]["selected_model"]]
    labels = ["Current full", "Selected"]
    metrics = ["auroc", "pr_auc", "brier"]
    colors = ["#64748b", "#14b8a6"]
    pieces = ['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="430" viewBox="0 0 900 430">',
              '<rect width="100%" height="100%" fill="#07111f"/>',
              '<text x="40" y="30" fill="#e2e8f0" font-family="sans-serif" font-size="18">Locked 2020+ evaluation</text>']
    for metric_index, metric in enumerate(metrics):
        x0 = 90 + metric_index * 270
        pieces.append(f'<text x="{x0+85}" y="385" text-anchor="middle" fill="#cbd5e1" font-family="sans-serif" font-size="14">{metric}</text>')
        for model_index, (name, label, color) in enumerate(zip(names, labels, colors)):
            value = float(results["evaluation"][name]["point"][metric])
            height = 280 * value
            x = x0 + model_index * 95
            y = 345 - height
            pieces.append(f'<rect x="{x}" y="{y}" width="70" height="{height}" rx="5" fill="{color}"/>')
            pieces.append(f'<text x="{x+35}" y="{y-7}" text-anchor="middle" fill="#e2e8f0" font-family="sans-serif" font-size="12">{value:.3f}</text>')
            if metric_index == 0:
                pieces.append(f'<text x="{x+35}" y="410" text-anchor="middle" fill="{color}" font-family="sans-serif" font-size="11">{label}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def _report(results: dict[str, Any]) -> str:
    selected = results["development_selection"]["selected_model"]
    before = results["evaluation"]["current_full_logistic"]
    after = results["evaluation"][selected]
    np_result = results["natural_products"]["metrics"]
    np_current = results["natural_products"]["current_full_metrics"]
    cluster = results["study_bias"]["CHEMBL5729139"]
    cluster_current = results["study_bias"]["CHEMBL5729139_current_full"]
    after_slope = after["point"]["calibration_slope"]
    after_slope_text = "not identifiable" if after_slope is None else f"{after_slope:.3f}"
    return f"""# Study-grouped assay-disagreement generalization analysis

## Decision protocol

The 2020+ evaluation was already summarized in the preceding experiment, so it is locked but not analyst-naive. No feature, candidate, hyperparameter, calibration mapping, or selection rule in this phase uses its outcomes. Selection uses only leave-one-study-out predictions from publications through {DEVELOPMENT_END_YEAR}, minimizing equal-study-weighted Brier score. Natural-product status is never a predictor.

The prespecified outcome remains a ≥{OUTCOME_THRESHOLD_LOG10:.1f}-log10 difference between earlier and later exact-IC50 medians from document- and assay-disjoint records in one compatible context. It is an assay disagreement signal, not drug efficacy, safety, biological failure, or clinical outcome.

## Development selection

There were {results['dataset']['development_events']} development compounds from {results['dataset']['development_studies']} held-out studies; {results['dataset']['development_singleton_studies']} studies contained one compound. The selected method was `{selected}` by {results['development_selection']['criterion']}. All candidate LOSO metrics are in `development_model_metrics.csv`.

## Locked evaluation

| Model | AUROC | 95% study-cluster CI | PR-AUC | 95% study-cluster CI | Brier | Calibration slope |
|---|---:|---:|---:|---:|---:|---:|
| Current full logistic | {before['point']['auroc']:.3f} | {before['clustered_ci95']['auroc']} | {before['point']['pr_auc']:.3f} | {before['clustered_ci95']['pr_auc']} | {before['point']['brier']:.3f} | {before['point']['calibration_slope']:.3f} |
| Development-selected `{selected}` | {after['point']['auroc']:.3f} | {after['clustered_ci95']['auroc']} | {after['point']['pr_auc']:.3f} | {after['clustered_ci95']['pr_auc']} | {after['point']['brier']:.3f} | {after_slope_text} |

The paired selected-minus-current differences and clustered intervals are recorded in `results.json`. The selected model's status is **{results['validation_decision']['status']}**: {results['validation_decision']['reason']}

## Study bias

Binary disagreement ICC was {results['study_bias']['evaluation_binary_icc']}. CHEMBL5729139 contained {cluster['molecules']} compounds and {cluster['disagreements']} disagreements (observed rate {cluster['observed_rate']:.3f}); mean selected risk was {cluster['mean_predicted_risk']:.3f}, while the current full model's mean risk was {cluster_current['mean_predicted_risk']:.3f}. The current full model's AUROC/PR-AUC after excluding this study were {cluster_current['auroc_without_study']}/{cluster_current['pr_auc_without_study']}. This study is a protocol/batch shift the model cannot identify from earlier compound evidence alone.

## Natural products

The evaluation subgroup contained {np_result['point']['n']} structure-verified natural products and {np_result['point']['positives']} disagreements. The prior full model had AUROC {np_current['point']['auroc']} and PR-AUC {np_current['point']['pr_auc']}; the development-selected fallback had AUROC {np_result['point']['auroc']} (95% study-cluster CI {np_result['clustered_ci95']['auroc']}), PR-AUC {np_result['point']['pr_auc']} (CI {np_result['clustered_ci95']['pr_auc']}), and Brier score {np_result['point']['brier']}. {results['natural_products']['error_summary']['outside_descriptor_domain']} compounds were outside the development descriptor domain. The dominant overall EGFR study contained no verified natural products, while NP compounds had materially longer and more heterogeneous evidence histories. Compound errors and descriptor shifts are in dedicated CSVs.

## Interpretation

{results['scientific_conclusion']}

## Limitations

- Only 212 development compounds and 50 development disagreements were available after strict duplicate filtering; 112/141 development studies were singletons.
- Publication ID is an imperfect proxy for laboratory/protocol independence, and structured substrate, construct, detection, incubation, and laboratory variables remain sparse.
- The locked evaluation was observed in the previous phase, so this is not a fully analyst-blinded replication.
- One EGFR publication dominates positive evaluation outcomes; molecule-level confidence intervals would be anticonservative.
- Natural-product estimates remain small and target-confounded. NPASS/COCONUT annotations are subgroup labels only.
- Cross-organism generalization is not estimable from these five human-target contexts.
"""


def run_generalization_experiment(
    specifications: Path,
    output_dir: Path,
    *,
    bootstrap_iterations: int = 1000,
    prediction_interval_iterations: int = 200,
    seed: int = 2026,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    frame, input_manifest = build_dataset(specifications)
    development = frame[frame.partition == "development"].copy().reset_index(drop=True)
    evaluation = frame[frame.partition == "evaluation"].copy().reset_index(drop=True)
    if set(development.molecule_id) & set(evaluation.molecule_id):
        raise AssertionError("molecule leakage")
    selection, oof, calibrators = development_model_selection(development, seed=seed)
    fitted, predictions = fit_and_predict_evaluation(
        development, evaluation, selection, calibrators, seed=seed
    )
    selected = selection["selected_model"]
    selected_calibrator = calibrators.get(selected)
    lower, upper = _bootstrap_selected_intervals(
        development, evaluation, selected, selected_calibrator,
        iterations=prediction_interval_iterations, seed=seed,
    )
    predictions["selected_probability_ci95_lower"] = lower
    predictions["selected_probability_ci95_upper"] = upper

    evaluation_results: dict[str, Any] = {}
    evaluation_names = ["constant_prevalence", *CANDIDATES]
    for index, name in enumerate(evaluation_names):
        column = f"probability_{name}"
        evaluation_results[name] = {
            "point": _metric_values(predictions.disagreement_label, predictions[column]),
            "clustered_ci95": clustered_metric_intervals(
                predictions, column, iterations=bootstrap_iterations, seed=seed + index
            ),
        }
    selected_column = f"probability_{selected}"
    comparison = clustered_model_comparison(
        predictions, selected_column, "probability_current_full_logistic",
        iterations=bootstrap_iterations, seed=seed + 100,
    )
    study_rows = _study_table(predictions, selected_column)
    current_study_rows = _study_table(predictions, "probability_current_full_logistic")
    development_study_rows = _study_table(oof.rename(columns={f"probability_{selected}": "selected_probability"}), "selected_probability")
    cluster_row = next(row for row in study_rows if row["study"] == "CHEMBL5729139")
    current_cluster_row = next(
        row for row in current_study_rows if row["study"] == "CHEMBL5729139"
    )

    natural_errors, feature_shifts, natural_summary = _natural_product_analysis(development, predictions)
    natural_frame = predictions[predictions.is_verified_natural_product]
    natural_metrics = {
        "point": _metric_values(natural_frame.disagreement_label, natural_frame[selected_column]),
        "clustered_ci95": clustered_metric_intervals(
            natural_frame, selected_column, iterations=bootstrap_iterations, seed=seed + 200
        ),
    }
    natural_current_metrics = {
        "point": _metric_values(
            natural_frame.disagreement_label, natural_frame.probability_current_full_logistic
        ),
        "clustered_ci95": clustered_metric_intervals(
            natural_frame, "probability_current_full_logistic",
            iterations=bootstrap_iterations, seed=seed + 201,
        ),
    }
    selected_ci = evaluation_results[selected]["clustered_ci95"]
    selected_point = evaluation_results[selected]["point"]
    supported = bool(
        selected_ci.get("auroc") and selected_ci["auroc"][0] > 0.5
        and selected_ci.get("pr_auc") and selected_ci["pr_auc"][0] > float(evaluation.disagreement_label.mean())
        and selected_point["calibration_slope"] is not None
        and 0.5 <= float(selected_point["calibration_slope"]) <= 1.5
    )
    validation_decision = {
        "status": "supported_for_research_risk_estimation" if supported else "not_supported_for_reliable_prediction",
        "rule": "study-clustered AUROC CI > 0.5, PR-AUC CI > evaluation prevalence, and calibration slope 0.5-1.5",
        "reason": (
            "All prespecified discrimination and calibration conditions were met."
            if supported else "At least one prespecified study-generalization or calibration condition failed."
        ),
    }
    conclusion = (
        "Study-grouped development selection produced a model meeting the limited retrospective support rule, but external replication remains necessary."
        if supported else
        "Study-grouped selection did not establish reliable generalization. Study identity and protocol shift remain stronger determinants of disagreement than transferable compound-history signals."
    )
    results = {
        "generated_at": datetime.now(UTC).isoformat(),
        "seed": seed,
        "dataset": {
            "eligible_events": len(frame),
            "development_events": len(development),
            "development_disagreements": int(development.disagreement_label.sum()),
            "development_studies": int(development.heldout_study_key.nunique()),
            "development_singleton_studies": int((development.groupby("heldout_study_key").size() == 1).sum()),
            "evaluation_events": len(evaluation),
            "evaluation_disagreements": int(evaluation.disagreement_label.sum()),
            "evaluation_studies": int(evaluation.heldout_study_key.nunique()),
            "evaluation_natural_products": int(evaluation.is_verified_natural_product.sum()),
        },
        "outcome_definition": f"absolute earlier-vs-later exact IC50 median shift >= {OUTCOME_THRESHOLD_LOG10} log10",
        "development_selection": selection,
        "evaluation": evaluation_results,
        "selected_vs_current_full": comparison,
        "validation_decision": validation_decision,
        "study_bias": {
            "development_binary_icc": _binary_icc(development),
            "evaluation_binary_icc": _binary_icc(evaluation),
            "CHEMBL5729139": cluster_row,
            "CHEMBL5729139_current_full": current_cluster_row,
        },
        "natural_products": {
            "metrics": natural_metrics,
            "current_full_metrics": natural_current_metrics,
            "error_summary": natural_summary,
        },
        "feature_contract": {
            "enhanced_molecular": ENHANCED_MOLECULAR,
            "compact_history": COMPACT_HISTORY,
            "compact_metadata": COMPACT_METADATA,
            "prohibited": ["heldout values", "heldout assay/document identifiers", "future assay outcomes", "natural-product status"],
        },
        "input_manifest": input_manifest,
        "scientific_conclusion": conclusion,
    }

    _write_csv(output_dir / "dataset.csv", frame)
    _write_csv(output_dir / "development_loso_predictions.csv", oof)
    _write_csv(output_dir / "evaluation_predictions.csv", predictions)
    _write_csv(output_dir / "development_study_sensitivity.csv", development_study_rows)
    _write_csv(output_dir / "evaluation_study_sensitivity.csv", study_rows)
    _write_csv(output_dir / "evaluation_study_sensitivity_current_full.csv", current_study_rows)
    _write_csv(output_dir / "natural_product_errors.csv", natural_errors)
    _write_csv(output_dir / "natural_product_feature_shifts.csv", feature_shifts)
    metric_rows = []
    for name, values in selection["candidate_metrics"].items():
        metric_rows.append({"model": name, **{f"pooled_{k}": v for k, v in values["pooled"].items()}, **{f"study_balanced_{k}": v for k, v in values["study_balanced"].items()}})
    _write_csv(output_dir / "development_model_metrics.csv", metric_rows)
    artifact = {
        "selected_model": selected,
        "selection_criterion": selection["criterion"],
        "base_model": fitted[selected] if selected != "constant_prevalence" else None,
        "constant_prevalence": fitted["constant_prevalence"],
        "validation_status": validation_decision["status"],
        "feature_contract": results["feature_contract"],
        "limitations": "Research assay-disagreement artifact; not efficacy, safety, biological failure, or clinical prediction.",
    }
    joblib.dump(artifact, output_dir / "study_generalized_model.joblib")
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    _write_plot(output_dir / "comparison.svg", results)
    manifest = {
        "specifications": str(specifications.resolve()),
        "specifications_sha256": _sha256(specifications),
        "inputs": input_manifest["inputs"],
        "leakage_contract": "LOSO development selection; cross-fitted study-group calibration; locked 2020+ evaluation never used for selection",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    command = (
        f'python -m app.assay_disagreement_generalization --specs "{specifications.resolve()}" '
        f'--output-dir "{output_dir.resolve()}" --bootstrap-iterations {bootstrap_iterations} '
        f'--prediction-interval-iterations {prediction_interval_iterations} --seed {seed}'
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Study-grouped assay-disagreement generalization experiment")
    parser.add_argument("--specs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--prediction-interval-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    result = run_generalization_experiment(
        args.specs, args.output_dir, bootstrap_iterations=args.bootstrap_iterations,
        prediction_interval_iterations=args.prediction_interval_iterations, seed=args.seed,
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "selected_model": result["development_selection"]["selected_model"],
        "validation_decision": result["validation_decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
