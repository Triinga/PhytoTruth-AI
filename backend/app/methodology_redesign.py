from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import HuberRegressor, LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, mean_absolute_error, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .assay_disagreement import ALL_NUMERIC_FEATURES, MOLECULAR_FEATURES
from .assay_disagreement_generalization import COMPACT_HISTORY, COMPACT_METADATA, ENHANCED_MOLECULAR


SEED = 2026
THRESHOLD = 1.0
HISTORY_FEATURES = list(dict.fromkeys(COMPACT_HISTORY + COMPACT_METADATA + [
    "prior_measurements", "prior_iqr", "prior_range", "prior_standard_deviation",
    "prior_document_threshold_disagreement",
]))
PROTOCOL_FEATURES = [
    "prior_description_reported_fraction", "prior_parameter_reported_fraction",
    "prior_direct_field_completeness_mean", "prior_organism_reported",
    "prior_tissue_reported", "prior_cell_type_reported", "prior_parameters_reported",
    "prior_variant_reported", "prior_methodology_reported", "prior_detection_reported",
    "prior_biological_system_text_reported", "prior_endpoint_reported",
    "prior_temperature_reported", "prior_ph_reported", "prior_duration_reported",
]
FULL_FEATURES = list(dict.fromkeys(HISTORY_FEATURES + ENHANCED_MOLECULAR + PROTOCOL_FEATURES))
FORBIDDEN_TOKENS = ("heldout_", "comparison", "disagreement_label", "absolute_pactivity_shift", "signed_pactivity_shift")
PROBABILITY_METHODS = ("constant", "historical_beta", "legacy_logistic", "rf_history", "rf_protocol")
MAGNITUDE_METHODS = ("constant", "historical_beta", "huber_magnitude", "rf_history", "rf_protocol")


def assert_feature_contract(features: Sequence[str]) -> None:
    offenders = [feature for feature in features if any(token in feature for token in FORBIDDEN_TOKENS)]
    if offenders:
        raise ValueError(f"Post-outcome or held-out features are forbidden: {offenders}")
    permitted_molecular = set(ENHANCED_MOLECULAR) | set(MOLECULAR_FEATURES)
    non_prior = [feature for feature in features if feature not in permitted_molecular and not feature.startswith("prior_")]
    if non_prior:
        raise ValueError(f"Features must be molecular descriptors or prior-only evidence: {non_prior}")


def load_development_dataset(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = pd.read_csv(path)
    total = len(frame)
    development = frame[(frame["partition"] == "development") & (frame["heldout_year"] <= 2019)].copy()
    development = development.dropna(subset=["absolute_pactivity_shift", "disagreement_label"])
    development["disagreement_label"] = development["disagreement_label"].astype(int)
    development["is_verified_natural_product"] = development["is_verified_natural_product"].astype(bool)
    if (development["absolute_pactivity_shift"] < 0).any():
        raise ValueError("Disagreement magnitudes must be non-negative")
    expected = (development["absolute_pactivity_shift"] >= THRESHOLD).astype(int)
    if not np.array_equal(expected.to_numpy(), development["disagreement_label"].to_numpy()):
        raise ValueError("Binary labels are inconsistent with the frozen 1.0 pChEMBL threshold")
    assert_feature_contract(FULL_FEATURES)
    assert_feature_contract(ALL_NUMERIC_FEATURES)
    missing = sorted(set(FULL_FEATURES + list(ALL_NUMERIC_FEATURES)) - set(development.columns))
    if missing:
        raise ValueError(f"Missing required development features: {missing}")
    manifest = {
        "source_rows": total,
        "development_rows": len(development),
        "excluded_locked_evaluation_rows": int(total - len(development)),
        "molecules": int(development.molecule_id.nunique()),
        "studies": int(development.heldout_study_key.nunique()),
        "targets": int(development.target_chembl_id.nunique()),
        "disagreements": int(development.disagreement_label.sum()),
        "natural_products": int(development.is_verified_natural_product.sum()),
        "natural_product_disagreements": int(development.loc[development.is_verified_natural_product, "disagreement_label"].sum()),
        "year_min": int(development.heldout_year.min()),
        "year_max": int(development.heldout_year.max()),
    }
    return development.reset_index(drop=True), manifest


def _connected_groups(frame: pd.DataFrame) -> np.ndarray:
    """Connected study/compound components prevent either entity crossing a fold."""
    parent: dict[str, str] = {}

    def find(value: str) -> str:
        parent.setdefault(value, value)
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: str, right: str) -> None:
        lroot, rroot = find(left), find(right)
        if lroot != rroot:
            parent[max(lroot, rroot)] = min(lroot, rroot)

    for row in frame.itertuples():
        union(f"study:{row.heldout_study_key}", f"molecule:{row.molecule_id}")
    return np.asarray([find(f"study:{study}") for study in frame.heldout_study_key], dtype=str)


def make_splits(frame: pd.DataFrame) -> dict[str, list[tuple[np.ndarray, np.ndarray]]]:
    components = _connected_groups(frame)
    grouped = list(GroupKFold(n_splits=5).split(frame, frame.disagreement_label, components))
    temporal: list[tuple[np.ndarray, np.ndarray]] = []
    for train_end, test_start, test_end in ((2010, 2011, 2013), (2013, 2014, 2016), (2016, 2017, 2019)):
        train = frame.index[frame.heldout_year <= train_end].to_numpy()
        test = frame.index[frame.heldout_year.between(test_start, test_end)].to_numpy()
        test_molecules = set(frame.loc[test, "molecule_id"])
        train = np.asarray([index for index in train if frame.loc[index, "molecule_id"] not in test_molecules], dtype=int)
        if len(train) and len(test):
            temporal.append((train, test))
    target: list[tuple[np.ndarray, np.ndarray]] = []
    for target_id in sorted(frame.target_chembl_id.unique()):
        test = frame.index[frame.target_chembl_id == target_id].to_numpy()
        test_molecules = set(frame.loc[test, "molecule_id"])
        train = frame.index[(frame.target_chembl_id != target_id) & ~frame.molecule_id.isin(test_molecules)].to_numpy()
        target.append((train, test))
    splits = {"study_grouped": grouped, "forward_temporal": temporal, "leave_target_out": target}
    validate_splits(frame, splits)
    return splits


def validate_splits(frame: pd.DataFrame, splits: dict[str, list[tuple[np.ndarray, np.ndarray]]]) -> None:
    for scheme, folds in splits.items():
        for train, test in folds:
            if set(train) & set(test):
                raise AssertionError(f"Row leakage in {scheme}")
            if set(frame.loc[train, "molecule_id"]) & set(frame.loc[test, "molecule_id"]):
                raise AssertionError(f"Compound leakage in {scheme}")
            if scheme == "study_grouped" and set(frame.loc[train, "heldout_study_key"]) & set(frame.loc[test, "heldout_study_key"]):
                raise AssertionError("Publication leakage in study-grouped validation")
            if scheme == "forward_temporal" and frame.loc[train, "heldout_year"].max() >= frame.loc[test, "heldout_year"].min():
                raise AssertionError("Temporal leakage")


def _preprocessor(features: Sequence[str], *, scale: bool) -> ColumnTransformer:
    steps: list[tuple[str, Any]] = [("impute", SimpleImputer(strategy="median"))]
    if scale:
        steps.append(("scale", StandardScaler()))
    return ColumnTransformer([
        ("numeric", Pipeline(steps), list(features)),
        ("target", OneHotEncoder(handle_unknown="ignore", sparse_output=False), ["target_chembl_id"]),
    ])


def _pipeline(estimator: Any, features: Sequence[str], *, scale: bool) -> Pipeline:
    return Pipeline([("preprocess", _preprocessor(features, scale=scale)), ("model", estimator)])


def _fit_fold(train: pd.DataFrame, test: pd.DataFrame, seed: int) -> dict[str, np.ndarray]:
    y_binary = train.disagreement_label.to_numpy(int)
    y_magnitude = train.absolute_pactivity_shift.to_numpy(float)
    prevalence = float(y_binary.mean())
    median_magnitude = float(np.median(y_magnitude))
    prior_documents = np.maximum(train.prior_documents.to_numpy(float), 0)
    outputs: dict[str, np.ndarray] = {
        "risk_constant": np.full(len(test), prevalence),
        "magnitude_constant": np.full(len(test), median_magnitude),
        "risk_historical_beta": (
            0.5 + test.prior_historical_tenfold_fraction.to_numpy(float) * np.maximum(test.prior_documents.to_numpy(float), 0)
        ) / (1.0 + np.maximum(test.prior_documents.to_numpy(float), 0)),
        "magnitude_historical_beta": np.maximum(test.prior_document_median_range.to_numpy(float), 0),
    }
    del prior_documents
    if len(np.unique(y_binary)) == 2:
        legacy = _pipeline(
            LogisticRegression(C=1.0, solver="liblinear", max_iter=2000, random_state=seed),
            ALL_NUMERIC_FEATURES, scale=True,
        )
        legacy.fit(train[list(ALL_NUMERIC_FEATURES) + ["target_chembl_id"]], y_binary)
        outputs["risk_legacy_logistic"] = legacy.predict_proba(test[list(ALL_NUMERIC_FEATURES) + ["target_chembl_id"]])[:, 1]
    else:
        outputs["risk_legacy_logistic"] = np.full(len(test), prevalence)
    huber = _pipeline(HuberRegressor(epsilon=1.35, alpha=0.0001, max_iter=10000), HISTORY_FEATURES, scale=True)
    huber.fit(train[HISTORY_FEATURES + ["target_chembl_id"]], np.log1p(y_magnitude))
    outputs["magnitude_huber_magnitude"] = np.clip(
        np.expm1(huber.predict(test[HISTORY_FEATURES + ["target_chembl_id"]])), 0, float(np.max(y_magnitude))
    )
    outputs["risk_huber_magnitude"] = outputs["magnitude_huber_magnitude"]
    for offset, (name, features) in enumerate((("rf_history", HISTORY_FEATURES), ("rf_protocol", FULL_FEATURES))):
        model = _pipeline(RandomForestRegressor(
            n_estimators=400, min_samples_leaf=5, max_features=0.7, bootstrap=True,
            random_state=seed + offset, n_jobs=-1,
        ), features, scale=False)
        model.fit(train[features + ["target_chembl_id"]], y_magnitude)
        transformed = model.named_steps["preprocess"].transform(test[features + ["target_chembl_id"]])
        tree_predictions = np.column_stack([tree.predict(transformed) for tree in model.named_steps["model"].estimators_])
        tree_predictions = np.maximum(tree_predictions, 0)
        outputs[f"magnitude_{name}"] = np.median(tree_predictions, axis=1)
        outputs[f"risk_{name}"] = np.mean(tree_predictions >= THRESHOLD, axis=1)
        outputs[f"lower_{name}"] = np.quantile(tree_predictions, 0.10, axis=1)
        outputs[f"upper_{name}"] = np.quantile(tree_predictions, 0.90, axis=1)
    return outputs


def cross_validated_predictions(frame: pd.DataFrame, splits: dict[str, list[tuple[np.ndarray, np.ndarray]]], seed: int) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    identity = [
        "molecule_id", "target_chembl_id", "heldout_study_key", "heldout_year",
        "disagreement_label", "absolute_pactivity_shift", "is_verified_natural_product",
        "canonical_smiles", "prior_document_ids", "prior_assay_ids",
        "heldout_document_ids", "heldout_assay_ids",
    ]
    for scheme, folds in splits.items():
        for fold, (train_indices, test_indices) in enumerate(folds):
            train, test = frame.iloc[train_indices], frame.iloc[test_indices]
            predicted = test[identity].copy()
            predicted.insert(0, "source_index", test.index)
            predicted.insert(0, "fold", fold)
            predicted.insert(0, "validation_scheme", scheme)
            for column, values in _fit_fold(train, test, seed + fold * 31 + len(rows) * 1009).items():
                predicted[column] = values
            rows.append(predicted)
    return pd.concat(rows, ignore_index=True)


def _risk_metrics(group: pd.DataFrame, column: str, *, probability: bool) -> dict[str, float | None]:
    y, score = group.disagreement_label.to_numpy(int), group[column].to_numpy(float)
    result: dict[str, float | None] = {
        "auroc": float(roc_auc_score(y, score)) if len(np.unique(y)) == 2 else None,
        "pr_auc": float(average_precision_score(y, score)) if len(np.unique(y)) == 2 else None,
        "brier": None, "ece_10": None,
    }
    if probability:
        probability_values = np.clip(score, 0, 1)
        result["brier"] = float(brier_score_loss(y, probability_values))
        bins = np.minimum((probability_values * 10).astype(int), 9)
        result["ece_10"] = float(sum(
            np.mean(bins == index) * abs(float(y[bins == index].mean()) - float(probability_values[bins == index].mean()))
            for index in range(10) if np.any(bins == index)
        ))
    return result


def _bootstrap_ci(group: pd.DataFrame, method: str, *, iterations: int, seed: int) -> dict[str, list[float] | None]:
    clusters = [cluster.index.to_numpy() for _, cluster in group.reset_index(drop=True).groupby("heldout_study_key", sort=True)]
    reset = group.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    sampled: dict[str, list[float]] = defaultdict(list)
    for _ in range(iterations):
        indices = np.concatenate([clusters[index] for index in rng.integers(0, len(clusters), len(clusters))])
        replicate = reset.iloc[indices]
        if f"risk_{method}" in replicate:
            risk = _risk_metrics(replicate, f"risk_{method}", probability=method in PROBABILITY_METHODS)
            for metric in ("auroc", "pr_auc", "brier", "ece_10"):
                if risk[metric] is not None and math.isfinite(float(risk[metric])):
                    sampled[metric].append(float(risk[metric]))
        if f"magnitude_{method}" in replicate:
            sampled["mae"].append(float(mean_absolute_error(replicate.absolute_pactivity_shift, replicate[f"magnitude_{method}"])))
    return {metric: ([float(value) for value in np.quantile(values, [0.025, 0.975])] if values else None) for metric, values in sampled.items()}


def evaluate(predictions: pd.DataFrame, *, iterations: int, seed: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    methods = sorted(set(PROBABILITY_METHODS + MAGNITUDE_METHODS + ("huber_magnitude",)))
    rows: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}
    for scheme, scheme_frame in predictions.groupby("validation_scheme", sort=True):
        nested[scheme] = {}
        for scope, scoped in (("all", scheme_frame), ("natural_products", scheme_frame[scheme_frame.is_verified_natural_product])):
            nested[scheme][scope] = {}
            for number, method in enumerate(methods):
                row: dict[str, Any] = {
                    "validation_scheme": scheme, "scope": scope, "method": method,
                    "n": len(scoped), "positives": int(scoped.disagreement_label.sum()),
                    "prevalence": float(scoped.disagreement_label.mean()) if len(scoped) else None,
                }
                if len(scoped) and f"risk_{method}" in scoped:
                    row.update(_risk_metrics(scoped, f"risk_{method}", probability=method in PROBABILITY_METHODS))
                else:
                    row.update({"auroc": None, "pr_auc": None, "brier": None, "ece_10": None})
                if len(scoped) and f"magnitude_{method}" in scoped:
                    row["mae"] = float(mean_absolute_error(scoped.absolute_pactivity_shift, scoped[f"magnitude_{method}"]))
                else:
                    row["mae"] = None
                if len(scoped) and f"lower_{method}" in scoped:
                    row["interval_80_coverage"] = float(np.mean(
                        (scoped.absolute_pactivity_shift >= scoped[f"lower_{method}"])
                        & (scoped.absolute_pactivity_shift <= scoped[f"upper_{method}"])
                    ))
                    row["interval_mean_width"] = float(np.mean(scoped[f"upper_{method}"] - scoped[f"lower_{method}"]))
                else:
                    row["interval_80_coverage"] = None
                    row["interval_mean_width"] = None
                intervals = _bootstrap_ci(scoped, method, iterations=iterations, seed=seed + number) if len(scoped) else {}
                for metric in ("auroc", "pr_auc", "brier", "ece_10", "mae"):
                    row[f"{metric}_ci95"] = json.dumps(intervals.get(metric))
                rows.append(row)
                nested[scheme][scope][method] = {key: value for key, value in row.items() if key not in {"validation_scheme", "scope", "method"}}
    return pd.DataFrame(rows), nested


def paired_protocol_ablation(predictions: pd.DataFrame, *, iterations: int, seed: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Paired cluster bootstrap: protocol+molecular RF minus history-only RF."""
    rows: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}
    for scheme, scheme_frame in predictions.groupby("validation_scheme", sort=True):
        nested[scheme] = {}
        for scope, scoped in (("all", scheme_frame), ("natural_products", scheme_frame[scheme_frame.is_verified_natural_product])):
            reset = scoped.reset_index(drop=True)
            clusters = [group.index.to_numpy() for _, group in reset.groupby("heldout_study_key", sort=True)]
            point_full = _risk_metrics(reset, "risk_rf_protocol", probability=True)
            point_history = _risk_metrics(reset, "risk_rf_history", probability=True)
            point = {
                "auroc": None if point_full["auroc"] is None or point_history["auroc"] is None else point_full["auroc"] - point_history["auroc"],
                "pr_auc": None if point_full["pr_auc"] is None or point_history["pr_auc"] is None else point_full["pr_auc"] - point_history["pr_auc"],
                "mae": float(mean_absolute_error(reset.absolute_pactivity_shift, reset.magnitude_rf_protocol) - mean_absolute_error(reset.absolute_pactivity_shift, reset.magnitude_rf_history)),
            }
            rng = np.random.default_rng(seed + len(rows))
            sampled: dict[str, list[float]] = defaultdict(list)
            for _ in range(iterations):
                indices = np.concatenate([clusters[index] for index in rng.integers(0, len(clusters), len(clusters))])
                replicate = reset.iloc[indices]
                full = _risk_metrics(replicate, "risk_rf_protocol", probability=True)
                history = _risk_metrics(replicate, "risk_rf_history", probability=True)
                for metric in ("auroc", "pr_auc"):
                    if full[metric] is not None and history[metric] is not None:
                        sampled[metric].append(float(full[metric] - history[metric]))
                sampled["mae"].append(float(
                    mean_absolute_error(replicate.absolute_pactivity_shift, replicate.magnitude_rf_protocol)
                    - mean_absolute_error(replicate.absolute_pactivity_shift, replicate.magnitude_rf_history)
                ))
            intervals = {metric: ([float(value) for value in np.quantile(values, [0.025, 0.975])] if values else None) for metric, values in sampled.items()}
            record = {"validation_scheme": scheme, "scope": scope}
            for metric in ("auroc", "pr_auc", "mae"):
                record[f"delta_{metric}"] = point[metric]
                record[f"delta_{metric}_ci95"] = json.dumps(intervals.get(metric))
            rows.append(record)
            nested[scheme][scope] = {"point": point, "study_clustered_ci95": intervals}
    return pd.DataFrame(rows), nested


def _fit_final_candidates(frame: pd.DataFrame, seed: int) -> dict[str, Any]:
    y_magnitude = frame.absolute_pactivity_shift.to_numpy(float)
    y_binary = frame.disagreement_label.to_numpy(int)
    fitted: dict[str, Any] = {
        "constant": {"prevalence": float(y_binary.mean()), "median_magnitude": float(np.median(y_magnitude))},
        "feature_contract": {"history": HISTORY_FEATURES, "protocol": FULL_FEATURES},
    }
    legacy = _pipeline(LogisticRegression(C=1.0, solver="liblinear", max_iter=2000, random_state=seed), ALL_NUMERIC_FEATURES, scale=True)
    legacy.fit(frame[list(ALL_NUMERIC_FEATURES) + ["target_chembl_id"]], y_binary)
    fitted["legacy_logistic"] = legacy
    huber = _pipeline(HuberRegressor(epsilon=1.35, alpha=0.0001, max_iter=10000), HISTORY_FEATURES, scale=True)
    huber.fit(frame[HISTORY_FEATURES + ["target_chembl_id"]], np.log1p(y_magnitude))
    fitted["huber_magnitude"] = {"model": huber, "prediction_bounds": [0.0, float(np.max(y_magnitude))]}
    for offset, (name, features) in enumerate((("rf_history", HISTORY_FEATURES), ("rf_protocol", FULL_FEATURES))):
        model = _pipeline(RandomForestRegressor(
            n_estimators=400, min_samples_leaf=5, max_features=0.7, bootstrap=True,
            random_state=seed + offset, n_jobs=-1,
        ), features, scale=False)
        model.fit(frame[features + ["target_chembl_id"]], y_magnitude)
        fitted[name] = model
    return fitted


def _plot(metrics: pd.DataFrame, path: Path) -> None:
    subset = metrics[(metrics.validation_scheme == "study_grouped") & (metrics.scope == "all")]
    methods = ["constant", "historical_beta", "legacy_logistic", "huber_magnitude", "rf_history", "rf_protocol"]
    colors = ["#64748b", "#0ea5e9", "#f59e0b", "#a78bfa", "#14b8a6", "#22c55e"]
    pieces = ['<svg xmlns="http://www.w3.org/2000/svg" width="1040" height="520" viewBox="0 0 1040 520">', '<rect width="1040" height="520" fill="#07111f"/>', '<text x="55" y="42" fill="#f8fafc" font-family="sans-serif" font-size="22">Development-only study-grouped validation</text>']
    for panel, metric in enumerate(("auroc", "mae")):
        x0 = 55 + panel * 500
        pieces.append(f'<text x="{x0}" y="78" fill="#cbd5e1" font-family="sans-serif" font-size="15">{metric.upper()}</text>')
        for index, (method, color) in enumerate(zip(methods, colors)):
            row = subset[subset.method == method]
            value = None if row.empty else row.iloc[0][metric]
            if pd.isna(value):
                continue
            y = 105 + index * 58
            width = max(0, min(390, float(value) * (390 if metric == "auroc" else 260)))
            pieces.append(f'<text x="{x0}" y="{y+17}" fill="#cbd5e1" font-family="sans-serif" font-size="12">{method.replace("_", " ")}</text>')
            pieces.append(f'<rect x="{x0+145}" y="{y}" width="{width}" height="24" rx="4" fill="{color}"/>')
            pieces.append(f'<text x="{x0+150+width}" y="{y+17}" fill="#f8fafc" font-family="sans-serif" font-size="12">{float(value):.3f}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def diagnostic_tables(frame: pd.DataFrame, predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    feature_rows: list[dict[str, Any]] = []
    for scope, scoped in (("all", frame), ("natural_products", frame[frame.is_verified_natural_product])):
        for feature in FULL_FEATURES:
            numeric = pd.to_numeric(scoped[feature], errors="coerce")
            feature_rows.append({
                "scope": scope, "feature": feature, "n": len(scoped),
                "missing_fraction": float(numeric.isna().mean()),
                "spearman_with_magnitude": float(numeric.corr(scoped.absolute_pactivity_shift, method="spearman")) if numeric.nunique(dropna=True) > 1 else None,
                "positive_median": float(numeric[scoped.disagreement_label == 1].median()) if (scoped.disagreement_label == 1).any() else None,
                "negative_median": float(numeric[scoped.disagreement_label == 0].median()) if (scoped.disagreement_label == 0).any() else None,
            })
    study_predictions = predictions[predictions.validation_scheme == "study_grouped"].copy()
    study_predictions["rf_protocol_absolute_error"] = np.abs(
        study_predictions.absolute_pactivity_shift - study_predictions.magnitude_rf_protocol
    )
    natural_errors = study_predictions[study_predictions.is_verified_natural_product].copy()
    natural_errors["error_type"] = np.select(
        [
            (natural_errors.disagreement_label == 1) & (natural_errors.risk_rf_protocol < 0.5),
            (natural_errors.disagreement_label == 0) & (natural_errors.risk_rf_protocol >= 0.5),
        ], ["missed_disagreement", "false_alarm"], default="correct_at_0.5",
    )
    natural_errors = natural_errors.sort_values(
        ["error_type", "rf_protocol_absolute_error", "molecule_id"], ascending=[True, False, True]
    )
    study_counts = frame.groupby("heldout_study_key", sort=True).agg(
        events=("molecule_id", "size"), disagreements=("disagreement_label", "sum"),
        natural_products=("is_verified_natural_product", "sum"),
    ).reset_index()
    diagnostics = {
        "singleton_studies": int((study_counts.events == 1).sum()),
        "studies": len(study_counts),
        "largest_study_events": int(study_counts.events.max()),
        "largest_study_disagreements": int(study_counts.loc[study_counts.events.idxmax(), "disagreements"]),
        "top_five_studies_fraction_of_events": float(study_counts.nlargest(5, "events").events.sum() / len(frame)),
        "natural_product_missed_disagreements_at_0_5": int((natural_errors.error_type == "missed_disagreement").sum()),
        "natural_product_false_alarms_at_0_5": int((natural_errors.error_type == "false_alarm").sum()),
    }
    return pd.DataFrame(feature_rows), natural_errors, diagnostics


def _report(results: dict[str, Any]) -> str:
    data = results["dataset"]
    metrics = results["metrics"]
    lines = [
        "# Assay Disagreement Methodology Redesign", "", "## Scope and leakage controls", "",
        f"This development-only analysis used {data['development_rows']} pre-2020 events ({data['disagreements']} disagreements) from {data['studies']} publication groups and {data['targets']} targets. It excluded {data['excluded_locked_evaluation_rows']} rows in the locked 2020+ partition. The separately assembled 90-case and 94-case independent datasets were not loaded.", "",
        "The primary outcome is absolute held-out pChEMBL shift. Binary ten-fold disagreement is secondary. Every validation fold separates molecules; study-grouped folds also separate publication groups, forward folds train strictly before their test years, and leave-target-out folds hold out a complete target.", "",
        "## Why the previous method failed", "",
        f"The original development set has {results['diagnostics']['singleton_studies']} singleton publication groups out of {results['diagnostics']['studies']}; consequently, study-specific shifts cannot be learned robustly. The earlier investigation also measured study-level label clustering (ICC 0.288) and found that its calibrated candidates were worse than a constant-prevalence model in development. Protocol enrichment did not recover the central missing conditions: across 1,050 assays, cell type and strain were never structured, while temperature, pH, methodology, endpoint, and parameter types were reported or detected in only 2.8%, 4.2%, 12.6%, 11.7%, and 1.2% respectively. Reporting indicators therefore cannot substitute for actual protocol values.", "",
        "## Development validation", "",
        "| Scheme | Scope | Method | n | AUROC | PR-AUC | MAE | Brier |", "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for scheme in ("study_grouped", "forward_temporal", "leave_target_out"):
        for scope in ("all", "natural_products"):
            for method in ("constant", "historical_beta", "legacy_logistic", "huber_magnitude", "rf_history", "rf_protocol"):
                row = metrics[scheme][scope][method]
                fmt = lambda value: "—" if value is None or pd.isna(value) else f"{float(value):.3f}"
                lines.append(f"| {scheme} | {scope} | {method} | {row['n']} | {fmt(row['auroc'])} | {fmt(row['pr_auc'])} | {fmt(row['mae'])} | {fmt(row['brier'])} |")
    decision = results["readiness"]
    ablation = results["paired_protocol_ablation"]["study_grouped"]
    all_delta = ablation["all"]
    natural_delta = ablation["natural_products"]
    lines += [
        "", "## Scientific interpretation", "",
        results["conclusion"], "",
        f"In the paired study-cluster bootstrap, adding molecular and prior protocol-completeness features changed overall AUROC by {all_delta['point']['auroc']:.3f} (95% CI {all_delta['study_clustered_ci95']['auroc']}) and MAE by {all_delta['point']['mae']:.3f} (95% CI {all_delta['study_clustered_ci95']['mae']}; negative favors the richer model). For natural products, the AUROC change was {natural_delta['point']['auroc']:.3f} (95% CI {natural_delta['study_clustered_ci95']['auroc']}). These intervals determine whether the added feature group has defensible value; point estimates alone do not.", "",
        f"**Readiness decision: `{decision['status']}`.** {decision['reason']}", "",
        "Random-forest tree quantiles are descriptive ensemble intervals, not validated confidence or prediction intervals. The historical-beta score is a Jeffreys-smoothed prior contradiction frequency, not a validated probability. Natural-product status was never used as a model feature.", "",
        f"At the fixed 0.5 descriptive RF vote threshold, the study-grouped natural-product analysis missed {results['diagnostics']['natural_product_missed_disagreements_at_0_5']} disagreements and produced {results['diagnostics']['natural_product_false_alarms_at_0_5']} false alarms. Compound-level provenance and residuals are in `natural_product_error_analysis.csv`; this threshold was not tuned.", "",
        "## Limitations", "",
        "- The development set is small, with only 50 disagreements and 14 natural-product disagreements; subgroup confidence intervals are consequently wide.",
        "- Protocol features mainly encode whether metadata were reported, not the missing experimental conditions themselves.",
        "- Leave-target-out validation is severe with only five targets, while forward-temporal folds cover only the later part of the development era.",
        "- Disagreement indicates a large difference between comparable database measurements; it does not establish which experiment is correct or imply clinical failure.",
        "- No result here updates or reinterprets any previous independent evaluation.",
    ]
    return "\n".join(lines) + "\n"


def run_experiment(dataset: Path, output_dir: Path, *, bootstrap_iterations: int = 1000, seed: int = SEED) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    frame, manifest = load_development_dataset(dataset)
    splits = make_splits(frame)
    predictions = cross_validated_predictions(frame, splits, seed)
    metric_table, nested_metrics = evaluate(predictions, iterations=bootstrap_iterations, seed=seed)
    ablation_table, ablation = paired_protocol_ablation(predictions, iterations=bootstrap_iterations, seed=seed + 7000)
    feature_table, natural_errors, diagnostics = diagnostic_tables(frame, predictions)
    primary = nested_metrics["study_grouped"]["all"]["rf_protocol"]
    temporal = nested_metrics["forward_temporal"]["all"]["rf_protocol"]
    history_np = nested_metrics["study_grouped"]["natural_products"]["rf_history"]
    protocol_np = nested_metrics["study_grouped"]["natural_products"]["rf_protocol"]
    criteria = {
        "study_grouped_auroc_above_chance": primary["auroc"] is not None and primary["auroc"] > 0.5,
        "study_grouped_pr_above_prevalence": primary["pr_auc"] is not None and primary["pr_auc"] > primary["prevalence"],
        "study_grouped_mae_below_constant": primary["mae"] < nested_metrics["study_grouped"]["all"]["constant"]["mae"],
        "forward_temporal_auroc_above_chance": temporal["auroc"] is not None and temporal["auroc"] > 0.5,
        "natural_product_not_degraded": protocol_np["auroc"] is not None and history_np["auroc"] is not None and protocol_np["auroc"] >= history_np["auroc"],
    }
    ready = all(criteria.values())
    readiness = {
        "status": "ready_for_new_independent_test" if ready else "not_ready_for_new_independent_test",
        "criteria": criteria,
        "reason": "All prespecified development criteria were met." if ready else "At least one prespecified generalization criterion failed; additional independent test use is not scientifically justified yet.",
    }
    conclusion = (
        "Magnitude modeling was evaluated because it retains more information than binary thresholding. "
        f"In study-grouped validation, the protocol RF had AUROC {primary['auroc']:.3f}, PR-AUC {primary['pr_auc']:.3f}, and MAE {primary['mae']:.3f}; "
        f"the constant baseline MAE was {nested_metrics['study_grouped']['all']['constant']['mae']:.3f}. "
        f"For natural products, protocol-RF AUROC was {protocol_np['auroc']} versus {history_np['auroc']} for history-only RF. "
        "These are development estimates, not independent validation or evidence of deployment readiness."
    )
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "protocol_version": "methodology-redesign-v1", "seed": seed,
        "dataset": manifest,
        "validation": {name: {"folds": len(folds), "test_predictions": int(sum(len(test) for _, test in folds))} for name, folds in splits.items()},
        "metrics": nested_metrics, "paired_protocol_ablation": ablation, "diagnostics": diagnostics,
        "readiness": readiness, "conclusion": conclusion,
        "independent_evaluation_datasets_loaded": False,
    }
    predictions.to_csv(output_dir / "development_predictions.csv", index=False)
    metric_table.to_csv(output_dir / "development_metrics.csv", index=False)
    ablation_table.to_csv(output_dir / "protocol_ablation.csv", index=False)
    feature_table.to_csv(output_dir / "feature_diagnostics.csv", index=False)
    natural_errors.to_csv(output_dir / "natural_product_error_analysis.csv", index=False)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_report(results), encoding="utf-8")
    _plot(metric_table, output_dir / "development_comparison.svg")
    joblib.dump({
        "artifact_version": "methodology-redesign-v1", "trained_on": "development partition through 2019 only",
        "not_independently_validated": True, "threshold": THRESHOLD, "models": _fit_final_candidates(frame, seed),
    }, output_dir / "candidate_models.joblib")
    manifest_out = {
        "files": sorted(path.name for path in output_dir.iterdir() if path.is_file()),
        "development_source": str(dataset), "rows": len(frame), "locked_rows_excluded": manifest["excluded_locked_evaluation_rows"],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest_out, indent=2), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Development-only assay disagreement methodology redesign")
    parser.add_argument("--dataset", type=Path, default=Path("experiments/assay_protocol_enrichment/enriched_assay_disagreement_events.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("experiments/assay_disagreement_methodology_redesign"))
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    results = run_experiment(args.dataset, args.output_dir, bootstrap_iterations=args.bootstrap_iterations, seed=args.seed)
    print(json.dumps({"dataset": results["dataset"], "readiness": results["readiness"]}, indent=2))


if __name__ == "__main__":
    main()
