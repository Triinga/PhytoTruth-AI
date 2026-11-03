from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Sequence

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.special import expit, logit
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .features import DESCRIPTOR_NAMES, bemis_murcko_scaffold, molecular_descriptors
from .multi_target_benchmark import load_target_specs
from .reversal_experiment import _jsonl, _sha256, natural_product_membership, select_compatible_measurements


OUTCOME_THRESHOLD_LOG10 = 1.0
ACTIVE_THRESHOLD = 6.0
DEVELOPMENT_END_YEAR = 2019
EVALUATION_START_YEAR = 2020

MOLECULAR_FEATURES = list(DESCRIPTOR_NAMES)
HISTORY_FEATURES = [
    "prior_median_pactivity", "prior_distance_to_active_threshold", "prior_measurements",
    "prior_documents", "prior_assays", "prior_year_span", "prior_mad", "prior_iqr",
    "prior_range", "prior_standard_deviation", "prior_document_median_range",
    "prior_document_threshold_disagreement",
]
PROVENANCE_FEATURES = [
    "prior_records_per_document", "prior_assays_per_document", "prior_journals",
    "prior_descriptions", "prior_missing_description_fraction", "prior_missing_journal_fraction",
]
ALL_NUMERIC_FEATURES = MOLECULAR_FEATURES + HISTORY_FEATURES + PROVENANCE_FEATURES

MODEL_FEATURES = {
    "constant_prevalence": [],
    "target_only": [],
    "molecular_only": MOLECULAR_FEATURES,
    "prior_potency": ["prior_median_pactivity", "prior_distance_to_active_threshold"],
    "evidence_history": HISTORY_FEATURES,
    "assay_provenance": PROVENANCE_FEATURES,
    "full": ALL_NUMERIC_FEATURES,
    "full_without_molecular": HISTORY_FEATURES + PROVENANCE_FEATURES,
    "full_without_history": MOLECULAR_FEATURES + PROVENANCE_FEATURES,
    "full_without_provenance": MOLECULAR_FEATURES + HISTORY_FEATURES,
}


def _mad(values: Sequence[float]) -> float:
    center = median(values)
    return float(median(abs(value - center) for value in values))


def _enhanced_molecular_descriptors(smiles: str) -> dict[str, float]:
    """Small interpretable descriptor extension; no learned or assay-outcome features."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES")
    heavy_atoms = max(1, molecule.GetNumHeavyAtoms())
    aromatic_atoms = sum(atom.GetIsAromatic() for atom in molecule.GetAtoms())
    hetero_atoms = sum(atom.GetAtomicNum() not in {1, 6} for atom in molecule.GetAtoms())
    return {
        "aromatic_atom_fraction": float(aromatic_atoms / heavy_atoms),
        "heteroatom_fraction": float(hetero_atoms / heavy_atoms),
        "bertz_complexity": float(Descriptors.BertzCT(molecule)),
    }


def _metric_values(
    labels: Sequence[int], probabilities: Sequence[float], *, calibration_diagnostics: bool = True
) -> dict[str, float | None]:
    y = np.asarray(labels, dtype=int)
    p = np.clip(np.asarray(probabilities, dtype=float), 1e-6, 1 - 1e-6)
    result: dict[str, float | None] = {
        "n": int(len(y)),
        "positives": int(y.sum()),
        "prevalence": float(y.mean()) if len(y) else None,
        "auroc": None,
        "pr_auc": None,
        "brier": None,
        "log_loss": None,
        "ece_10": None,
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    if not len(y):
        return result
    result["brier"] = float(brier_score_loss(y, p))
    result["log_loss"] = float(log_loss(y, p, labels=[0, 1]))
    bins = np.minimum((p * 10).astype(int), 9)
    result["ece_10"] = float(sum(
        np.mean(bins == index) * abs(float(y[bins == index].mean()) - float(p[bins == index].mean()))
        for index in range(10) if np.any(bins == index)
    ))
    if len(np.unique(y)) == 2:
        result["auroc"] = float(roc_auc_score(y, p))
        result["pr_auc"] = float(average_precision_score(y, p))
        if calibration_diagnostics and float(np.ptp(p)) > 1e-12:
            calibration = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
            calibration.fit(logit(p).reshape(-1, 1), y)
            result["calibration_intercept"] = float(calibration.intercept_[0])
            result["calibration_slope"] = float(calibration.coef_[0, 0])
    return result


def clustered_metric_intervals(
    frame: pd.DataFrame,
    probability_column: str,
    *,
    iterations: int,
    seed: int,
) -> dict[str, list[float] | None]:
    """Resample held-out study clusters within target; never individual activity rows."""
    cluster_indices = {
        str(target): [cluster.index.to_numpy() for _, cluster in group.groupby("heldout_study_key", sort=True)]
        for target, group in frame.reset_index(drop=True).groupby("target_chembl_id", sort=True)
    }
    frame = frame.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    sampled_metrics: dict[str, list[float]] = defaultdict(list)
    for _ in range(iterations):
        sampled_indices = []
        for target in sorted(cluster_indices):
            clusters = cluster_indices[target]
            choices = rng.integers(0, len(clusters), size=len(clusters))
            sampled_indices.extend(clusters[int(choice)] for choice in choices)
        replicate = frame.iloc[np.concatenate(sampled_indices)]
        metrics = _metric_values(
            replicate["disagreement_label"], replicate[probability_column], calibration_diagnostics=False
        )
        for name in ("auroc", "pr_auc", "brier", "log_loss", "ece_10"):
            if metrics[name] is not None and math.isfinite(float(metrics[name])):
                sampled_metrics[name].append(float(metrics[name]))
    return {
        name: ([float(v) for v in np.quantile(values, [0.025, 0.975])] if values else None)
        for name, values in sampled_metrics.items()
    }


def clustered_model_comparison(
    frame: pd.DataFrame,
    candidate_column: str,
    reference_column: str,
    *,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    """Paired study-cluster bootstrap for model metric differences."""
    point_candidate = _metric_values(frame.disagreement_label, frame[candidate_column])
    point_reference = _metric_values(frame.disagreement_label, frame[reference_column])
    metric_names = ("auroc", "pr_auc", "brier")
    point = {
        name: (float(point_candidate[name]) - float(point_reference[name]))
        if point_candidate[name] is not None and point_reference[name] is not None else None
        for name in metric_names
    }
    frame = frame.reset_index(drop=True)
    cluster_indices = {
        str(target): [cluster.index.to_numpy() for _, cluster in group.groupby("heldout_study_key", sort=True)]
        for target, group in frame.groupby("target_chembl_id", sort=True)
    }
    rng = np.random.default_rng(seed)
    differences: dict[str, list[float]] = defaultdict(list)
    for _ in range(iterations):
        sampled_indices = []
        for target in sorted(cluster_indices):
            clusters = cluster_indices[target]
            choices = rng.integers(0, len(clusters), size=len(clusters))
            sampled_indices.extend(clusters[int(choice)] for choice in choices)
        sample = frame.iloc[np.concatenate(sampled_indices)]
        candidate = _metric_values(sample.disagreement_label, sample[candidate_column], calibration_diagnostics=False)
        reference = _metric_values(sample.disagreement_label, sample[reference_column], calibration_diagnostics=False)
        for name in metric_names:
            if candidate[name] is not None and reference[name] is not None:
                differences[name].append(float(candidate[name]) - float(reference[name]))
    return {
        name: {
            "difference": point[name],
            "clustered_ci95": [float(value) for value in np.quantile(differences[name], [0.025, 0.975])]
            if differences[name] else None,
        }
        for name in metric_names
    }


def _case_from_records(
    target: dict[str, Any],
    molecule_id: str,
    records: Sequence[dict[str, Any]],
    natural: dict[str, list[dict[str, Any]]],
) -> tuple[dict[str, Any] | None, str | None]:
    years = sorted({int((row.get("raw_source") or {}).get("document_year")) for row in records})
    if len(years) < 2:
        return None, "single_publication_year"
    heldout_year = years[-1]
    prior = [row for row in records if int(row["raw_source"]["document_year"]) < heldout_year]
    heldout = [row for row in records if int(row["raw_source"]["document_year"]) == heldout_year]
    prior_docs = {str(row.get("document_chembl_id")) for row in prior}
    heldout_docs = {str(row.get("document_chembl_id")) for row in heldout}
    prior_assays = {str(row.get("assay_chembl_id")) for row in prior}
    heldout_assays = {str(row.get("assay_chembl_id")) for row in heldout}
    if prior_docs & heldout_docs:
        return None, "document_overlap"
    if prior_assays & heldout_assays:
        return None, "assay_overlap"
    smiles_values = {str(row.get("canonical_smiles")) for row in records if row.get("canonical_smiles")}
    if len(smiles_values) != 1:
        return None, "missing_or_inconsistent_structure"
    smiles = next(iter(smiles_values))
    try:
        descriptors = molecular_descriptors(smiles)
        scaffold = bemis_murcko_scaffold(smiles)
    except ValueError:
        return None, "invalid_structure"
    prior_values = [float(row["computed_pactivity"]) for row in prior]
    heldout_values = [float(row["computed_pactivity"]) for row in heldout]
    document_values: dict[str, list[float]] = defaultdict(list)
    for row in prior:
        document_values[str(row["document_chembl_id"])].append(float(row["computed_pactivity"]))
    doc_medians = [median(values) for values in document_values.values()]
    assay_values: dict[str, list[float]] = defaultdict(list)
    year_values: dict[int, list[float]] = defaultdict(list)
    for row in prior:
        assay_values[str(row["assay_chembl_id"])].append(float(row["computed_pactivity"]))
        year_values[int(row["raw_source"]["document_year"])].append(float(row["computed_pactivity"]))
    assay_medians = [median(values) for values in assay_values.values()]
    ordered_years = sorted(year_values)
    year_medians = [median(year_values[year]) for year in ordered_years]
    prior_center = float(median(prior_values))
    heldout_center = float(median(heldout_values))
    delta = abs(heldout_center - prior_center)
    descriptions = [(row.get("raw_source") or {}).get("assay_description") for row in prior]
    journals = [(row.get("raw_source") or {}).get("document_journal") for row in prior]
    provenance = natural.get(molecule_id, [])
    sources = sorted({str(item.get("source")) for item in provenance if item.get("source")})
    case: dict[str, Any] = {
        "target_chembl_id": target["target_chembl_id"],
        "target_name": target["target_name"],
        "molecule_id": molecule_id,
        "canonical_smiles": smiles,
        "murcko_scaffold": scaffold,
        "heldout_year": heldout_year,
        "partition": "development" if heldout_year <= DEVELOPMENT_END_YEAR else "evaluation",
        "heldout_study_key": "+".join(sorted(heldout_docs)),
        "prior_activity_ids": ";".join(sorted(str(row["activity_id"]) for row in prior)),
        "heldout_activity_ids": ";".join(sorted(str(row["activity_id"]) for row in heldout)),
        "prior_document_ids": ";".join(sorted(prior_docs)),
        "heldout_document_ids": ";".join(sorted(heldout_docs)),
        "prior_assay_ids": ";".join(sorted(prior_assays)),
        "heldout_assay_ids": ";".join(sorted(heldout_assays)),
        "prior_median_pactivity": prior_center,
        "heldout_median_pactivity": heldout_center,
        "signed_pactivity_shift": heldout_center - prior_center,
        "absolute_pactivity_shift": delta,
        "disagreement_label": int(delta >= OUTCOME_THRESHOLD_LOG10),
        "active_threshold_contradiction": int((prior_center >= ACTIVE_THRESHOLD) != (heldout_center >= ACTIVE_THRESHOLD)),
        "is_verified_natural_product": bool(provenance),
        "natural_product_sources": ";".join(sources),
        "natural_product_provenance_count": len(provenance),
        "prior_distance_to_active_threshold": abs(prior_center - ACTIVE_THRESHOLD),
        "prior_measurements": len(prior),
        "prior_documents": len(prior_docs),
        "prior_assays": len(prior_assays),
        "prior_year_span": max(int(row["raw_source"]["document_year"]) for row in prior) - min(int(row["raw_source"]["document_year"]) for row in prior),
        "prior_mad": _mad(prior_values),
        "prior_iqr": float(np.quantile(prior_values, 0.75) - np.quantile(prior_values, 0.25)),
        "prior_range": max(prior_values) - min(prior_values),
        "prior_standard_deviation": float(np.std(prior_values)),
        "prior_document_median_range": max(doc_medians) - min(doc_medians),
        "prior_document_threshold_disagreement": float(np.mean([(value >= ACTIVE_THRESHOLD) != (prior_center >= ACTIVE_THRESHOLD) for value in doc_medians])),
        "prior_records_per_document": len(prior) / len(prior_docs),
        "prior_assays_per_document": len(prior_assays) / len(prior_docs),
        "prior_journals": len({str(value) for value in journals if value}),
        "prior_descriptions": len({str(value) for value in descriptions if value}),
        "prior_missing_description_fraction": sum(not value for value in descriptions) / len(descriptions),
        "prior_missing_journal_fraction": sum(not value for value in journals) / len(journals),
    }
    case.update(descriptors)
    case.update(_enhanced_molecular_descriptors(smiles))
    case.update({
        "prior_assay_median_range": max(assay_medians) - min(assay_medians),
        "prior_historical_tenfold_fraction": float(np.mean([
            abs(value - prior_center) >= OUTCOME_THRESHOLD_LOG10 for value in doc_medians
        ])),
        "prior_recent_shift": float(year_medians[-1] - median(year_medians[:-1])) if len(year_medians) > 1 else 0.0,
        "prior_temporal_slope": float(np.polyfit(ordered_years, year_medians, 1)[0]) if len(ordered_years) > 1 else 0.0,
        "prior_cell_keyword_fraction": float(np.mean([
            "cell" in str(value).lower() for value in descriptions if value
        ])) if any(descriptions) else 0.0,
        "prior_recombinant_keyword_fraction": float(np.mean([
            "recombinant" in str(value).lower() for value in descriptions if value
        ])) if any(descriptions) else 0.0,
        "prior_substrate_keyword_fraction": float(np.mean([
            "substrate" in str(value).lower() for value in descriptions if value
        ])) if any(descriptions) else 0.0,
    })
    return case, None


def build_dataset(specifications: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    target_summaries: list[dict[str, Any]] = []
    input_manifest: list[dict[str, Any]] = []
    total_exclusions: Counter[str] = Counter()
    for spec in load_target_specs(specifications):
        activities = Path(spec["activities"])
        natural_path = Path(spec["natural_products"])
        all_records = list(_jsonl(activities))
        quality_exclusions: Counter[str] = Counter()
        quality_records = []
        for record in all_records:
            if record.get("data_validity_comment"):
                quality_exclusions["data_validity_comment"] += 1
            elif record.get("potential_duplicate"):
                quality_exclusions["potential_duplicate"] += 1
            elif record.get("pchembl_consistent") is False:
                quality_exclusions["pchembl_normalization_mismatch"] += 1
            else:
                quality_records.append(record)
        selected, context_key, context_exclusions = select_compatible_measurements(quality_records)
        natural = natural_product_membership(natural_path)
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in selected:
            grouped[str(row["molecule_chembl_id"])].append(row)
        target_exclusions: Counter[str] = Counter()
        before = len(cases)
        for molecule_id in sorted(grouped):
            case, reason = _case_from_records(spec, molecule_id, grouped[molecule_id], natural)
            if reason:
                target_exclusions[reason] += 1
            else:
                cases.append(case or {})
        total_exclusions.update(target_exclusions)
        target_cases = cases[before:]
        target_summaries.append({
            "target_chembl_id": spec["target_chembl_id"],
            "target_name": spec["target_name"],
            "selected_context_key": context_key,
            "selected_measurements": len(selected),
            "selected_molecules": len(grouped),
            "eligible_events": len(target_cases),
            "development_events": sum(row["partition"] == "development" for row in target_cases),
            "evaluation_events": sum(row["partition"] == "evaluation" for row in target_cases),
            "evaluation_disagreements": sum(row["partition"] == "evaluation" and row["disagreement_label"] for row in target_cases),
            "evaluation_natural_products": sum(row["partition"] == "evaluation" and row["is_verified_natural_product"] for row in target_cases),
            "context_exclusions": dict(context_exclusions),
            "quality_exclusions": dict(quality_exclusions),
            "event_exclusions": dict(target_exclusions),
        })
        total_exclusions.update(quality_exclusions)
        input_manifest.append({
            "target_chembl_id": spec["target_chembl_id"],
            "activities": str(activities), "activities_sha256": _sha256(activities),
            "natural_products": str(natural_path), "natural_products_sha256": _sha256(natural_path),
        })
    frame = pd.DataFrame(cases).sort_values(["target_chembl_id", "heldout_year", "molecule_id"]).reset_index(drop=True)
    if frame.empty:
        raise ValueError("No eligible independent later-study events")
    development = set(frame.loc[frame.partition == "development", "molecule_id"])
    evaluation = set(frame.loc[frame.partition == "evaluation", "molecule_id"])
    overlapping_molecules = development & evaluation
    if overlapping_molecules:
        removed = int(((frame.partition == "development") & frame.molecule_id.isin(overlapping_molecules)).sum())
        frame = frame[~((frame.partition == "development") & frame.molecule_id.isin(overlapping_molecules))].reset_index(drop=True)
        total_exclusions["cross_target_molecule_overlap_removed_from_development"] += removed
    development = set(frame.loc[frame.partition == "development", "molecule_id"])
    evaluation = set(frame.loc[frame.partition == "evaluation", "molecule_id"])
    if development & evaluation:
        raise AssertionError("molecule leakage between development and evaluation")
    if not (frame.loc[frame.partition == "development", "heldout_year"] <= DEVELOPMENT_END_YEAR).all():
        raise AssertionError("development year leakage")
    if not (frame.loc[frame.partition == "evaluation", "heldout_year"] >= EVALUATION_START_YEAR).all():
        raise AssertionError("evaluation year leakage")
    for summary in target_summaries:
        target_frame = frame[frame.target_chembl_id == summary["target_chembl_id"]]
        summary["eligible_events_before_cross_target_deduplication"] = summary.pop("eligible_events")
        summary["eligible_events"] = len(target_frame)
        summary["development_events"] = int((target_frame.partition == "development").sum())
        summary["evaluation_events"] = int((target_frame.partition == "evaluation").sum())
        summary["evaluation_disagreements"] = int(target_frame.loc[target_frame.partition == "evaluation", "disagreement_label"].sum())
        summary["evaluation_natural_products"] = int(target_frame.loc[target_frame.partition == "evaluation", "is_verified_natural_product"].sum())
    return frame, {
        "inputs": input_manifest,
        "targets": target_summaries,
        "exclusions": dict(total_exclusions),
        "outcome": f"absolute difference between prior-study and latest-year median exact IC50 pActivity >= {OUTCOME_THRESHOLD_LOG10} log10 unit",
        "split": f"development held-out publications <= {DEVELOPMENT_END_YEAR}; untouched evaluation >= {EVALUATION_START_YEAR}",
    }


def _pipeline(numeric_features: Sequence[str]) -> Pipeline:
    transformers = []
    if numeric_features:
        transformers.append(("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ]), list(numeric_features)))
    transformers.append(("target", OneHotEncoder(handle_unknown="ignore"), ["target_chembl_id"]))
    return Pipeline([
        ("preprocess", ColumnTransformer(transformers)),
        ("classifier", LogisticRegression(C=1.0, solver="liblinear", max_iter=2000, random_state=2026)),
    ])


def fit_models(development: pd.DataFrame, evaluation: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    predictions = evaluation.copy()
    models: dict[str, Any] = {}
    prevalence = float(development["disagreement_label"].mean())
    predictions["probability_constant_prevalence"] = prevalence
    models["constant_prevalence"] = {"prevalence": prevalence}
    for name, features in MODEL_FEATURES.items():
        if name == "constant_prevalence":
            continue
        model = _pipeline(features)
        model.fit(development[list(features) + ["target_chembl_id"]], development["disagreement_label"])
        predictions[f"probability_{name}"] = model.predict_proba(evaluation[list(features) + ["target_chembl_id"]])[:, 1]
        models[name] = model
    return models, predictions


def bootstrap_prediction_intervals(
    development: pd.DataFrame,
    evaluation: pd.DataFrame,
    *,
    iterations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    estimates: list[np.ndarray] = []
    for iteration in range(iterations):
        pieces = []
        for _, group in development.groupby("target_chembl_id", sort=True):
            indices = rng.integers(0, len(group), len(group))
            pieces.append(group.iloc[indices])
        sample = pd.concat(pieces, ignore_index=True)
        if sample["disagreement_label"].nunique() < 2:
            continue
        model = _pipeline(MODEL_FEATURES["full"])
        columns = MODEL_FEATURES["full"] + ["target_chembl_id"]
        model.fit(sample[columns], sample["disagreement_label"])
        estimates.append(model.predict_proba(evaluation[columns])[:, 1])
    matrix = np.asarray(estimates)
    return np.quantile(matrix, 0.025, axis=0), np.quantile(matrix, 0.975, axis=0)


def _model_metrics(predictions: pd.DataFrame, *, iterations: int, seed: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}
    for index, name in enumerate(MODEL_FEATURES):
        probability = f"probability_{name}"
        point = _metric_values(predictions.disagreement_label, predictions[probability])
        intervals = clustered_metric_intervals(predictions, probability, iterations=iterations, seed=seed + index)
        nested[name] = {"all": {"point": point, "clustered_ci95": intervals}, "targets": {}, "natural_products": {}}
        rows.append({"model": name, "subgroup": "all", **point, **{f"{key}_ci95": json.dumps(value) for key, value in intervals.items()}})
        for target, group in predictions.groupby("target_chembl_id", sort=True):
            metrics = _metric_values(group.disagreement_label, group[probability])
            nested[name]["targets"][str(target)] = metrics
            rows.append({"model": name, "subgroup": str(target), **metrics})
        natural = predictions[predictions.is_verified_natural_product]
        natural_metrics = _metric_values(natural.disagreement_label, natural[probability])
        nested[name]["natural_products"] = natural_metrics
        rows.append({"model": name, "subgroup": "verified_natural_products", **natural_metrics})
        novel = predictions[predictions.novel_scaffold]
        novel_metrics = _metric_values(novel.disagreement_label, novel[probability])
        nested[name]["novel_scaffolds"] = novel_metrics
        rows.append({"model": name, "subgroup": "novel_scaffolds", **novel_metrics})
    return rows, nested


def _coefficients(model: Pipeline) -> list[dict[str, Any]]:
    names = model.named_steps["preprocess"].get_feature_names_out()
    values = model.named_steps["classifier"].coef_[0]
    return sorted(
        ({"feature": str(name), "coefficient": float(value), "odds_ratio": float(math.exp(value))} for name, value in zip(names, values)),
        key=lambda row: abs(row["coefficient"]), reverse=True,
    )


def _write_csv(path: Path, rows: Iterable[dict[str, Any]] | pd.DataFrame) -> None:
    if isinstance(rows, pd.DataFrame):
        rows.to_csv(path, index=False, lineterminator="\n")
        return
    values = list(rows)
    if not values:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def _comparison_svg(path: Path, metrics: dict[str, Any]) -> None:
    models = ["constant_prevalence", "target_only", "molecular_only", "evidence_history", "assay_provenance", "full"]
    colors = ["#64748b", "#94a3b8", "#38bdf8", "#f59e0b", "#a78bfa", "#14b8a6"]
    values = [metrics[name]["all"]["point"]["pr_auc"] or 0 for name in models]
    width, height, left, top, chart_h = 900, 430, 80, 45, 300
    bar_w = 90
    gap = 45
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
              '<rect width="100%" height="100%" fill="#07111f"/>',
              '<text x="40" y="28" fill="#e2e8f0" font-family="sans-serif" font-size="18">Untouched 2020+ evaluation: PR-AUC</text>']
    for tick in range(6):
        value = tick / 5
        y = top + chart_h * (1 - value)
        pieces.append(f'<line x1="{left}" x2="860" y1="{y}" y2="{y}" stroke="#1e293b"/>')
        pieces.append(f'<text x="42" y="{y+5}" fill="#94a3b8" font-family="sans-serif" font-size="12">{value:.1f}</text>')
    for index, (name, value, color) in enumerate(zip(models, values, colors)):
        x = left + 30 + index * (bar_w + gap)
        bar_h = chart_h * value
        y = top + chart_h - bar_h
        pieces.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{bar_h}" rx="5" fill="{color}"/>')
        pieces.append(f'<text x="{x + bar_w/2}" y="{y-7}" text-anchor="middle" fill="#e2e8f0" font-family="sans-serif" font-size="12">{value:.3f}</text>')
        label = name.replace("_", " ")
        pieces.append(f'<text x="{x + bar_w/2}" y="{top+chart_h+22}" text-anchor="middle" fill="#cbd5e1" font-family="sans-serif" font-size="11">{label}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def _report(results: dict[str, Any]) -> str:
    data = results["dataset"]
    metrics = results["evaluation_metrics"]
    full = metrics["full"]["all"]
    baseline = metrics["constant_prevalence"]["all"]
    natural = metrics["full"]["natural_products"]
    target_rows = "\n".join(
        f"| {row['target_chembl_id']} | {row['development_events']} | {row['evaluation_events']} | {row['evaluation_disagreements']} | {row['evaluation_natural_products']} |"
        for row in results["manifest"]["targets"]
    )
    support = results["validation_decision"]
    return f"""# Evidence Transferability and Assay Disagreement Prediction

## Question and outcome

Can information already available for a compound predict whether an **independent, later ChEMBL study** will differ materially from its earlier exact-IC50 evidence? The primary disagreement outcome is an absolute shift of at least {OUTCOME_THRESHOLD_LOG10:.1f} log10 pActivity (ten-fold potency). This is an assay-reproducibility signal, not evidence of biological failure, efficacy, safety, or clinical outcome.

Only exact, normalized IC50 records in each target's prespecified compatible human single-protein context were used. Earlier and held-out activity, document, and assay identifiers are disjoint. Units, endpoint, target organism, assay type, and BAO format are therefore not pooled across incompatible contexts.

## Design

- Development outcomes: latest independent publication through {DEVELOPMENT_END_YEAR}; n={data['development_events']}, disagreements={data['development_disagreements']}.
- Untouched evaluation outcomes: latest independent publication from {EVALUATION_START_YEAR}; n={data['evaluation_events']}, disagreements={data['evaluation_disagreements']}.
- Molecules are disjoint across development and evaluation. Features use only structure and evidence preceding that molecule's held-out publication. Held-out assay metadata is never a predictor.
- The model is a prespecified L2 logistic regression. Comparators and ablations were fixed before evaluation. No evaluation-set tuning was performed.
- Intervals resample held-out study clusters within target; molecules from the same publication are not treated as independent studies.

| Target | Development | Evaluation | Evaluation disagreements | Evaluation NPs |
|---|---:|---:|---:|---:|
{target_rows}

## Independent evaluation

The constant development-prevalence baseline obtained AUROC {baseline['point']['auroc']} and PR-AUC {baseline['point']['pr_auc']:.3f}. The full model obtained AUROC {full['point']['auroc']:.3f} (study-clustered 95% CI {full['clustered_ci95']['auroc']}) and PR-AUC {full['point']['pr_auc']:.3f} (95% CI {full['clustered_ci95']['pr_auc']}). Its Brier score was {full['point']['brier']:.3f}, calibration intercept {full['point']['calibration_intercept']:.3f}, and calibration slope {full['point']['calibration_slope']:.3f}.

The verified-natural-product subgroup had n={natural['n']} and {natural['positives']} disagreements; full-model AUROC={natural['auroc']} and PR-AUC={natural['pr_auc']}. This subgroup is exploratory because study and positive-event counts are small.

**Validation decision:** `{support['status']}`. The prespecified support rule required the study-clustered AUROC interval to exclude 0.5 and the PR-AUC interval to exceed evaluation prevalence. {support['reason']}

Full numerical baselines, ablations, target-specific results, calibration, and uncertainty are in `model_metrics.csv` and `results.json`. `predictions.csv` preserves event provenance and per-compound bootstrap risk intervals. Coefficients are standardized associations, not causal effects.

## Transferability findings

{results['scientific_conclusion']}

The main source of uncertainty is study-level protocol shift: {results['study_cluster_analysis']['largest_cluster']} is the largest held-out publication cluster ({results['study_cluster_analysis']['largest_cluster_molecules']} compounds, {results['study_cluster_analysis']['largest_cluster_disagreements']} disagreements; {results['study_cluster_analysis']['fraction_of_all_disagreements_in_largest_cluster']:.1%} of all evaluation disagreements). A model that cannot see a future protocol cannot fully anticipate such a batch shift. The sensitivity outcome crossing pActivity {ACTIVE_THRESHOLD:.1f} is reported descriptively and is not substituted for the prespecified ten-fold outcome.

## Limitations

- ChEMBL assay descriptions do not provide uniformly structured protocol variables such as substrate concentration, incubation time, construct, detection technology, or laboratory; residual protocol heterogeneity remains.
- The split is publication-temporal but uses the current ChEMBL and NPASS/COCONUT releases. Natural-product status is used only for subgroup evaluation, never as a predictor.
- Latest-publication events are retrospective. Database deposition and duplicated literature values can still create selection effects despite document/assay disjointness and duplicate flags.
- Only five human targets and one compatible exact-IC50 context per target were studied. Cross-organism transferability is not estimable.
- Natural-product estimates are underpowered and cannot establish separate performance.
- Disagreement does not identify which study is correct and must not be interpreted as clinical failure.
"""


def run_experiment(
    specifications: Path,
    output_dir: Path,
    *,
    bootstrap_iterations: int = 1000,
    prediction_interval_iterations: int = 200,
    seed: int = 2026,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    frame, manifest = build_dataset(specifications)
    development = frame[frame.partition == "development"].copy()
    evaluation = frame[frame.partition == "evaluation"].copy()
    if development.disagreement_label.sum() < 20 or evaluation.disagreement_label.sum() < 20:
        raise ValueError("Insufficient positive disagreement events for supported modeling")
    development_scaffolds = set(development.murcko_scaffold)
    evaluation["novel_scaffold"] = ~evaluation.murcko_scaffold.isin(development_scaffolds)
    models, predictions = fit_models(development, evaluation)
    lower, upper = bootstrap_prediction_intervals(
        development, evaluation, iterations=prediction_interval_iterations, seed=seed
    )
    predictions["full_probability_ci95_lower"] = lower
    predictions["full_probability_ci95_upper"] = upper
    metric_rows, metric_results = _model_metrics(predictions, iterations=bootstrap_iterations, seed=seed)
    comparisons: dict[str, Any] = {}
    comparison_rows: list[dict[str, Any]] = []
    for index, reference in enumerate((
        "constant_prevalence", "molecular_only", "full_without_molecular",
        "full_without_history", "full_without_provenance",
    )):
        value = clustered_model_comparison(
            predictions, "probability_full", f"probability_{reference}",
            iterations=bootstrap_iterations, seed=seed + 100 + index,
        )
        comparisons[f"full_vs_{reference}"] = value
        for metric, item in value.items():
            comparison_rows.append({
                "comparison": f"full_vs_{reference}", "metric": metric,
                "difference": item["difference"], "clustered_ci95": json.dumps(item["clustered_ci95"]),
            })
    coefficients = _coefficients(models["full"])
    cluster_counts = predictions.groupby(["target_chembl_id", "heldout_study_key"], sort=True).agg(
        molecules=("molecule_id", "size"), disagreements=("disagreement_label", "sum"),
        natural_products=("is_verified_natural_product", "sum"), year=("heldout_year", "first"),
        mean_absolute_shift=("absolute_pactivity_shift", "mean"),
    ).reset_index()
    largest = cluster_counts.sort_values(["molecules", "heldout_study_key"], ascending=[False, True]).iloc[0]
    full_metrics = metric_results["full"]["all"]["point"]
    target_only = metric_results["target_only"]["all"]["point"]
    molecular = metric_results["molecular_only"]["all"]["point"]
    evidence = metric_results["evidence_history"]["all"]["point"]
    full_ci = metric_results["full"]["all"]["clustered_ci95"]
    supported = bool(
        full_ci.get("auroc") and full_ci["auroc"][0] > 0.5
        and full_ci.get("pr_auc") and full_ci["pr_auc"][0] > float(evaluation.disagreement_label.mean())
    )
    validation_decision = {
        "status": "supported_for_research_risk_estimation" if supported else "not_supported_for_reliable_prediction",
        "rule": "study-clustered AUROC CI lower bound > 0.5 and PR-AUC CI lower bound > untouched evaluation prevalence",
        "reason": (
            "Both prespecified criteria were met; this remains retrospective research validation."
            if supported else
            "At least one criterion failed, so the serialized model is an auditable research artifact and must not be presented as a validated predictor."
        ),
    }
    conclusion = (
        "The prespecified full model did not establish reliable prediction beyond simple baselines. " if not supported else
        "The prespecified full model met the limited retrospective support rule. "
    ) + (
        f"Full PR-AUC={full_metrics['pr_auc']:.3f}; target-only={target_only['pr_auc']:.3f}; molecular-only={molecular['pr_auc']:.3f}; earlier-evidence-only={evidence['pr_auc']:.3f}. "
        "Performance is reported as predictive association only; no evidence component is claimed transferable unless its ablation improves independently evaluated performance with compatible uncertainty."
    )
    results = {
        "generated_at": datetime.now(UTC).isoformat(),
        "seed": seed,
        "dataset": {
            "eligible_events": len(frame),
            "development_events": len(development),
            "development_disagreements": int(development.disagreement_label.sum()),
            "evaluation_events": len(evaluation),
            "evaluation_disagreements": int(evaluation.disagreement_label.sum()),
            "evaluation_prevalence": float(evaluation.disagreement_label.mean()),
            "evaluation_natural_products": int(evaluation.is_verified_natural_product.sum()),
            "evaluation_natural_product_disagreements": int(predictions.loc[predictions.is_verified_natural_product, "disagreement_label"].sum()),
            "evaluation_novel_scaffolds": int(evaluation.novel_scaffold.sum()),
        },
        "outcome_definition": manifest["outcome"],
        "leakage_contract": manifest["split"] + "; molecule, document, and assay disjointness enforced; held-out metadata excluded from predictors",
        "feature_groups": MODEL_FEATURES,
        "evaluation_metrics": metric_results,
        "paired_model_comparisons": comparisons,
        "validation_decision": validation_decision,
        "study_cluster_analysis": {
            "clusters": len(cluster_counts),
            "largest_cluster": str(largest.heldout_study_key),
            "largest_cluster_molecules": int(largest.molecules),
            "largest_cluster_disagreements": int(largest.disagreements),
            "fraction_of_all_disagreements_in_largest_cluster": float(largest.disagreements / evaluation.disagreement_label.sum()),
        },
        "manifest": manifest,
        "scientific_conclusion": conclusion,
    }
    _write_csv(output_dir / "dataset.csv", frame)
    prediction_columns = [column for column in predictions.columns if not column.startswith("probability_")] + [column for column in predictions.columns if column.startswith("probability_")]
    _write_csv(output_dir / "predictions.csv", predictions[prediction_columns])
    _write_csv(output_dir / "model_metrics.csv", metric_rows)
    _write_csv(output_dir / "model_comparisons.csv", comparison_rows)
    _write_csv(output_dir / "coefficients.csv", coefficients)
    _write_csv(output_dir / "heldout_study_clusters.csv", cluster_counts)
    joblib.dump({
        "model": models["full"], "model_name": "interpretable_l2_logistic_regression",
        "features": MODEL_FEATURES["full"], "outcome_threshold_log10": OUTCOME_THRESHOLD_LOG10,
        "development_end_year": DEVELOPMENT_END_YEAR, "seed": seed,
        "validation_status": validation_decision["status"],
        "limitations": "Research artifact for assay-disagreement risk only; not drug efficacy, safety, or clinical performance.",
    }, output_dir / "assay_disagreement_model.joblib")
    _comparison_svg(output_dir / "comparison.svg", metric_results)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "report.md").write_text(_report(results), encoding="utf-8")
    manifest_output = {
        **manifest,
        "specifications": str(specifications.resolve()),
        "specifications_sha256": _sha256(specifications),
        "code_contract": "all predictors derive from structure or evidence strictly preceding each held-out independent publication",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest_output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    command = (
        f'python -m app.assay_disagreement --specs "{specifications.resolve()}" --output-dir "{output_dir.resolve()}" '
        f'--bootstrap-iterations {bootstrap_iterations} --prediction-interval-iterations {prediction_interval_iterations} --seed {seed}'
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict independent later-study assay disagreement")
    parser.add_argument("--specs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--prediction-interval-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    result = run_experiment(
        args.specs, args.output_dir, bootstrap_iterations=args.bootstrap_iterations,
        prediction_interval_iterations=args.prediction_interval_iterations, seed=args.seed,
    )
    print(json.dumps({"output_dir": str(args.output_dir), "dataset": result["dataset"], "conclusion": result["scientific_conclusion"]}, indent=2))


if __name__ == "__main__":
    main()
