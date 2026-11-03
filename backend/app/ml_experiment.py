from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import numpy as np
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

from .features import bemis_murcko_scaffold, feature_config, molecular_feature_vector
from .ml import (
    SEED,
    available_model_names,
    make_regressor,
    max_tanimoto_similarity,
    nearest_neighbor_threshold,
    pack_fingerprints,
    regression_metrics,
)
from .reversal_experiment import _jsonl, natural_product_membership, select_compatible_measurements
from .structures import standardize_smiles


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_dataset(
    activities: Path,
    natural_products: Path,
    *,
    context_key: str | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Aggregate exact compatible measurements to one target per standardized structure."""
    natural = natural_product_membership(natural_products)
    selected, selected_context, exclusions = select_compatible_measurements(_jsonl(activities), context_key=context_key)
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    structures: dict[str, Any] = {}
    rejected_structures = 0
    for record in selected:
        smiles = record.get("canonical_smiles")
        if not smiles:
            rejected_structures += 1
            continue
        try:
            structure = standardize_smiles(str(smiles))
        except (ValueError, RuntimeError):
            rejected_structures += 1
            continue
        grouped[structure.standard_inchi_key].append(record)
        structures[structure.standard_inchi_key] = structure

    rows: list[dict[str, Any]] = []
    for inchi_key in sorted(grouped):
        records = grouped[inchi_key]
        structure = structures[inchi_key]
        values = [float(record["computed_pactivity"]) for record in records]
        molecule_ids = sorted({str(record["molecule_chembl_id"]) for record in records})
        provenance = [item for molecule_id in molecule_ids for item in natural.get(molecule_id, [])]
        rows.append({
            "structure_id": inchi_key,
            "canonical_smiles": structure.canonical_smiles,
            "molecule_chembl_ids": ";".join(molecule_ids),
            "target_pactivity": float(median(values)),
            "measurement_count": len(values),
            "pactivity_min": min(values),
            "pactivity_max": max(values),
            "first_document_year": min(int(record["document_year"]) for record in records),
            "last_document_year": max(int(record["document_year"]) for record in records),
            "is_verified_natural_product": bool(provenance),
            "natural_product_sources": ";".join(sorted({str(item.get("source")) for item in provenance})),
            "scaffold": bemis_murcko_scaffold(structure.canonical_smiles),
        })
    if len(rows) < 50:
        raise ValueError("Fewer than 50 unique compatible structures are available for modeling")
    metadata = {
        "assay_context_key": selected_context,
        "assay_context": selected[0]["assay_context"],
        "selected_measurements": len(selected),
        "unique_standardized_structures": len(rows),
        "verified_natural_products": sum(row["is_verified_natural_product"] for row in rows),
        "rejected_structures": rejected_structures,
        "normalization_exclusions": exclusions,
    }
    return rows, metadata


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _conformal_quantile(residuals: np.ndarray, coverage: float = 0.90) -> float:
    if not len(residuals):
        raise ValueError("Calibration residuals are empty")
    probability = min(1.0, np.ceil((len(residuals) + 1) * coverage) / len(residuals))
    return float(np.quantile(residuals, probability, method="higher"))


def assign_scaffold_folds(scaffolds: Iterable[str], folds: int = 5) -> np.ndarray:
    groups = np.asarray(list(scaffolds))
    if folds < 2 or len(np.unique(groups)) < folds:
        raise ValueError("Requested folds require at least that many unique scaffolds")
    assignments = np.empty(len(groups), dtype=int)
    splitter = GroupKFold(n_splits=folds)
    dummy = np.zeros(len(groups))
    for fold, (_, test) in enumerate(splitter.split(dummy, dummy, groups), start=1):
        assignments[test] = fold
    return assignments


def _cluster_bootstrap_ci(
    truth: np.ndarray,
    prediction: np.ndarray,
    scaffolds: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, list[float]]:
    rng = np.random.default_rng(seed)
    unique = np.unique(scaffolds)
    indices = {group: np.flatnonzero(scaffolds == group) for group in unique}
    samples = {name: [] for name in ("mae", "rmse", "r2")}
    for _ in range(iterations):
        chosen = rng.choice(unique, size=len(unique), replace=True)
        selected = np.concatenate([indices[group] for group in chosen])
        metrics = regression_metrics(truth[selected], prediction[selected])
        for name in samples:
            samples[name].append(metrics[name])
    return {name: [float(value) for value in np.quantile(values, [0.025, 0.975])] for name, values in samples.items()}


def _paired_cluster_bootstrap_difference(
    truth: np.ndarray,
    reference_prediction: np.ndarray,
    candidate_prediction: np.ndarray,
    scaffolds: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    unique = np.unique(scaffolds)
    indices = {group: np.flatnonzero(scaffolds == group) for group in unique}
    reference_point = regression_metrics(truth, reference_prediction)
    candidate_point = regression_metrics(truth, candidate_prediction)
    samples = {name: [] for name in ("mae", "rmse", "r2")}
    for _ in range(iterations):
        selected = np.concatenate([indices[group] for group in rng.choice(unique, size=len(unique), replace=True)])
        reference = regression_metrics(truth[selected], reference_prediction[selected])
        candidate = regression_metrics(truth[selected], candidate_prediction[selected])
        for name in samples:
            samples[name].append(candidate[name] - reference[name])
    return {
        name: {
            "reference": reference_point[name],
            "candidate": candidate_point[name],
            "difference_candidate_minus_reference": candidate_point[name] - reference_point[name],
            "difference_ci95": [float(value) for value in np.quantile(samples[name], [0.025, 0.975])],
        }
        for name in samples
    }


def _spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2:
        return 0.0
    left_rank = np.argsort(np.argsort(left, kind="stable"), kind="stable").astype(float)
    right_rank = np.argsort(np.argsort(right, kind="stable"), kind="stable").astype(float)
    if np.std(left_rank) == 0 or np.std(right_rank) == 0:
        return 0.0
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def _evaluate_model(
    name: str,
    X: np.ndarray,
    y: np.ndarray,
    scaffolds: np.ndarray,
    packed: np.ndarray,
    natural_mask: np.ndarray,
    fold_assignments: np.ndarray,
    *,
    seed: int,
    bootstrap_iterations: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    prediction = np.empty(len(y), dtype=float)
    lower = np.empty(len(y), dtype=float)
    upper = np.empty(len(y), dtype=float)
    similarity = np.empty(len(y), dtype=float)
    inside = np.empty(len(y), dtype=bool)
    fold_diagnostics: list[dict[str, Any]] = []
    for fold in sorted(np.unique(fold_assignments)):
        outer_test = np.flatnonzero(fold_assignments == fold)
        outer_train = np.flatnonzero(fold_assignments != fold)
        splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=seed + int(fold))
        fit_local, calibration_local = next(splitter.split(X[outer_train], y[outer_train], scaffolds[outer_train]))
        fit = outer_train[fit_local]
        calibration = outer_train[calibration_local]
        model = make_regressor(name, seed=seed + int(fold))
        model.fit(X[fit], y[fit])
        calibration_prediction = model.predict(X[calibration])
        half_width = _conformal_quantile(np.abs(y[calibration] - calibration_prediction))
        fold_prediction = np.asarray(model.predict(X[outer_test]), dtype=float)
        threshold = nearest_neighbor_threshold(packed[fit])
        fold_similarity = max_tanimoto_similarity(packed[outer_test], packed[fit])
        inflation = np.maximum(1.0, threshold / np.maximum(fold_similarity, 1e-6))
        fold_width = half_width * inflation
        prediction[outer_test] = fold_prediction
        lower[outer_test] = fold_prediction - fold_width
        upper[outer_test] = fold_prediction + fold_width
        similarity[outer_test] = fold_similarity
        inside[outer_test] = fold_similarity >= threshold
        fold_values = regression_metrics(y[outer_test], fold_prediction)
        fold_diagnostics.append({
            "model": name,
            "fold": int(fold),
            "fit_molecules": len(fit),
            "calibration_molecules": len(calibration),
            "test_molecules": len(outer_test),
            "fit_scaffolds": len(set(scaffolds[fit])),
            "calibration_scaffolds": len(set(scaffolds[calibration])),
            "test_scaffolds": len(set(scaffolds[outer_test])),
            **fold_values,
            "conformal_half_width_90": half_width,
            "ad_similarity_threshold": threshold,
            "outside_domain": int(np.count_nonzero(~inside[outer_test])),
        })
    overall = regression_metrics(y, prediction)
    confidence = _cluster_bootstrap_ci(y, prediction, scaffolds, iterations=bootstrap_iterations, seed=seed)
    covered = (y >= lower) & (y <= upper)
    absolute_error = np.abs(y - prediction)
    in_domain_metrics = regression_metrics(y[inside], prediction[inside]) if np.any(inside) else None
    outside_metrics = regression_metrics(y[~inside], prediction[~inside]) if np.any(~inside) else None
    natural_metrics = regression_metrics(y[natural_mask], prediction[natural_mask]) if np.any(natural_mask) else None
    result = {
        "overall": overall,
        "scaffold_cluster_bootstrap_ci95": confidence,
        "natural_product_subgroup": {"n": int(np.count_nonzero(natural_mask)), "metrics": natural_metrics},
        "uncertainty": {
            "target_coverage": 0.90,
            "empirical_coverage": float(np.mean(covered)),
            "mean_interval_width": float(np.mean(upper - lower)),
            "median_interval_width": float(np.median(upper - lower)),
            "absolute_error_interval_width_spearman": _spearman(absolute_error, upper - lower),
            "method": "nested scaffold split-conformal absolute residuals; distance-inflated outside domain",
        },
        "applicability_domain": {
            "inside_count": int(np.count_nonzero(inside)),
            "outside_count": int(np.count_nonzero(~inside)),
            "inside_metrics": in_domain_metrics,
            "outside_metrics": outside_metrics,
            "mean_test_to_train_tanimoto": float(np.mean(similarity)),
            "method": "5th percentile of training nearest-neighbor Morgan Tanimoto, derived within each fold",
        },
    }
    predictions = [{
        "row_index": index,
        "model": name,
        "fold": int(fold_assignments[index]),
        "observed_pactivity": float(y[index]),
        "predicted_pactivity": float(prediction[index]),
        "prediction_interval_low": float(lower[index]),
        "prediction_interval_high": float(upper[index]),
        "interval_contains_observation": bool(covered[index]),
        "nearest_training_tanimoto": float(similarity[index]),
        "inside_applicability_domain": bool(inside[index]),
    } for index in range(len(y))]
    return result, predictions, fold_diagnostics


def _comparison_svg(path: Path, results: dict[str, Any]) -> None:
    models = list(results["models"])
    metrics = ("mae", "rmse", "r2")
    colors = ("#95a9b8", "#42d392", "#4aa3ff")
    width, height = 1000, 500
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#07131f"/>',
        '<text x="45" y="42" fill="#e8f5ef" font-family="Arial" font-size="24" font-weight="700">Scaffold-held-out bioactivity prediction</text>',
        '<text x="45" y="67" fill="#95a9b8" font-family="Arial" font-size="13">Outer-fold predictions; lower is better for MAE/RMSE, higher is better for R²</text>',
    ]
    for panel, metric in enumerate(metrics):
        x0 = 55 + panel * 310
        values = [float(results["models"][name]["overall"][metric]) for name in models]
        if metric == "r2":
            maximum = max(1.0, max(values))
            minimum = min(0.0, min(values))
        else:
            maximum, minimum = max(values) * 1.15, 0.0
        parts.append(f'<text x="{x0}" y="105" fill="#d6e6df" font-family="Arial" font-size="16">{metric.upper()}</text>')
        for index, (name, value) in enumerate(zip(models, values)):
            scaled = (value - minimum) / max(maximum - minimum, 1e-9)
            bar_height = max(2.0, scaled * 250)
            x = x0 + 15 + index * 82
            y = 385 - bar_height
            parts.append(f'<rect x="{x}" y="{y:.1f}" width="55" height="{bar_height:.1f}" rx="4" fill="{colors[index % len(colors)]}"/>')
            parts.append(f'<text x="{x + 27.5}" y="{y - 7:.1f}" text-anchor="middle" fill="#e8f5ef" font-family="Arial" font-size="11">{value:.3f}</text>')
            parts.append(f'<text transform="translate({x + 12},420) rotate(35)" fill="#95a9b8" font-family="Arial" font-size="11">{name}</text>')
        parts.append(f'<line x1="{x0}" y1="385" x2="{x0 + 265}" y2="385" stroke="#385064"/>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def _report(results: dict[str, Any]) -> str:
    best = results["selection"]["best_learned_model"]
    baseline = results["models"]["mean_baseline"]["overall"]
    chosen = results["models"][best]
    status = "outperformed" if results["selection"]["beats_mean_baseline_on_mae"] else "did not outperform"
    lines = [
        "# PhytoTruth Molecular Bioactivity Model",
        "",
        f"Generated: {results['generated_at']}",
        "",
        "## Outcome",
        "",
        f"The selected learned model is **{best}**. It {status} the mean baseline on scaffold-held-out MAE.",
        "",
        "This model predicts experimental pActivity for the defined ChEMBL assay context. It does not predict clinical efficacy, safety, or a probability that a compound is a drug.",
        "",
        "## Dataset",
        "",
        f"- Target/context: `{results['dataset']['assay_context_key']}` — {json.dumps(results['dataset']['assay_context'], sort_keys=True)}.",
        f"- {results['dataset']['selected_measurements']} exact compatible measurements aggregated to {results['dataset']['unique_standardized_structures']} standardized structures.",
        f"- {results['dataset']['verified_natural_products']} structures are verified natural products through the supplied structure crossmatch.",
        f"- {results['validation']['unique_scaffolds']} Bemis–Murcko scaffold groups; {results['validation']['folds']} outer folds.",
        "- Replicate measurements were aggregated by standard InChIKey before splitting, preventing the same standardized structure from crossing folds.",
        "",
        "## Scaffold-held-out results",
        "",
        "| Model | MAE | RMSE | R² | Natural-product MAE | 90% interval coverage | Outside AD |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, result in results["models"].items():
        metrics = result["overall"]
        natural = result["natural_product_subgroup"]["metrics"]
        lines.append(
            f"| {name} | {metrics['mae']:.3f} | {metrics['rmse']:.3f} | {metrics['r2']:.3f} | "
            f"{natural['mae']:.3f} | {result['uncertainty']['empirical_coverage']:.3f} | "
            f"{result['applicability_domain']['outside_count']} |"
        )
    ci = chosen["scaffold_cluster_bootstrap_ci95"]
    learned_comparison = results["model_comparisons"].get("xgboost_vs_random_forest")
    lines.extend([
        "",
        f"Selected-model MAE 95% scaffold-cluster bootstrap CI: [{ci['mae'][0]:.3f}, {ci['mae'][1]:.3f}].",
        "",
        "## Validation and uncertainty",
        "",
        "- Outer test folds contain scaffolds absent from their model-fitting and calibration partitions.",
        "- Each outer training set is split again by scaffold into fitting and uncertainty-calibration partitions.",
        "- The 90% interval is based on calibration absolute residuals and inflated when Tanimoto similarity falls below the fold-specific applicability threshold.",
        "- Applicability thresholds are the fifth percentile of nearest-neighbor Morgan similarity within training data. This is an operational heuristic, not proof that an in-domain prediction is correct.",
        f"- Selected-model empirical 90% interval coverage was {chosen['uncertainty']['empirical_coverage']:.3f}; mean interval width was {chosen['uncertainty']['mean_interval_width']:.3f} pActivity units.",
        "",
        "## Model selection",
        "",
        f"The learned model with the lowest outer-fold MAE was selected. Mean-baseline MAE was {baseline['mae']:.3f}; selected-model MAE was {chosen['overall']['mae']:.3f}. No hyperparameter search used outer-fold outcomes.",
        "",
        (
            "XGBoost-minus-Random-Forest MAE was "
            f"{learned_comparison['mae']['difference_candidate_minus_reference']:.3f} "
            f"(95% paired scaffold-cluster bootstrap CI [{learned_comparison['mae']['difference_ci95'][0]:.3f}, "
            f"{learned_comparison['mae']['difference_ci95'][1]:.3f}])."
            if learned_comparison else
            "XGBoost was unavailable, so no XGBoost-versus-Random-Forest comparison was made."
        ),
        "",
        f"The interval-width versus absolute-error Spearman diagnostic was {chosen['uncertainty']['absolute_error_interval_width_spearman']:.3f}. Near-zero association means the per-molecule width did not meaningfully rank errors, even though aggregate coverage was close to its target.",
        "",
        "## Limitations",
        "",
        "- Model selection and performance estimation use the same fixed cross-validation run; an untouched external target-specific test set is not available.",
        "- ChEMBL measurements are heterogeneous literature evidence despite strict metadata filtering and per-structure aggregation.",
        "- Scaffold splits are deliberately difficult but do not reproduce every real deployment distribution shift.",
        "- Prediction intervals are cross-validated residual diagnostics, not guaranteed prospective coverage, especially outside the applicability domain.",
        "- Natural-product subgroup metrics are observational; natural-product status is not used as a predictive feature.",
        "- The Evidence Reversal Experiment remains inconclusive/negative and is not evidence that its reliability heuristic is validated.",
        "",
        "## Reproduction",
        "",
        "Use `command.txt`, `manifest.json`, `environment.json`, and `feature_config.json` in this directory. The serialized model is `best_model.joblib`.",
        "",
    ])
    return "\n".join(lines)


def run_training_experiment(
    activities: Path,
    natural_products: Path,
    output_dir: Path,
    *,
    context_key: str | None = None,
    folds: int = 5,
    bootstrap_iterations: int = 1000,
    seed: int = SEED,
) -> dict[str, Any]:
    import joblib
    import pandas
    import rdkit
    import scipy
    import sklearn
    from rdkit import RDLogger

    if folds < 3:
        raise ValueError("At least three scaffold folds are required")
    RDLogger.DisableLog("rdApp.warning")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, dataset_metadata = prepare_dataset(activities, natural_products, context_key=context_key)
    config = feature_config()
    X = np.stack([molecular_feature_vector(row["canonical_smiles"]) for row in rows]).astype(np.float32)
    y = np.asarray([row["target_pactivity"] for row in rows], dtype=float)
    scaffolds = np.asarray([row["scaffold"] for row in rows])
    natural_mask = np.asarray([row["is_verified_natural_product"] for row in rows], dtype=bool)
    packed = pack_fingerprints(X, int(config["n_bits"]))
    unique_scaffolds = np.unique(scaffolds)
    if len(unique_scaffolds) < folds:
        raise ValueError("Insufficient unique scaffolds for requested folds")
    fold_assignments = assign_scaffold_folds(scaffolds, folds)
    for index, row in enumerate(rows):
        row["fold"] = int(fold_assignments[index])
    _write_csv(output_dir / "training_dataset.csv", rows)

    models: dict[str, Any] = {}
    fold_rows: list[dict[str, Any]] = []
    prediction_rows: list[dict[str, Any]] = []
    for model_index, name in enumerate(available_model_names()):
        result, predictions, diagnostics = _evaluate_model(
            name,
            X,
            y,
            scaffolds,
            packed,
            natural_mask,
            fold_assignments,
            seed=seed + 1000 * model_index,
            bootstrap_iterations=bootstrap_iterations,
        )
        models[name] = result
        fold_rows.extend(diagnostics)
        # Add dataset identity to each prediction without exposing held-out values to fitting.
        for prediction in predictions:
            row = rows[int(prediction.pop("row_index"))]
            prediction_rows.append({
                "structure_id": row["structure_id"],
                "molecule_chembl_ids": row["molecule_chembl_ids"],
                "is_verified_natural_product": row["is_verified_natural_product"],
                "scaffold": row["scaffold"],
                **prediction,
            })

    learned = [name for name in models if name != "mean_baseline"]
    best = min(learned, key=lambda name: (models[name]["overall"]["mae"], name))
    baseline_mae = models["mean_baseline"]["overall"]["mae"]
    best_mae = models[best]["overall"]["mae"]
    prediction_lookup = {
        (item["model"], item["structure_id"]): item["predicted_pactivity"]
        for item in prediction_rows
    }
    predictions_by_model = {
        name: np.asarray([prediction_lookup[(name, row["structure_id"])] for row in rows], dtype=float)
        for name in models
    }
    comparisons = {
        "xgboost_vs_random_forest": _paired_cluster_bootstrap_difference(
            y,
            predictions_by_model["random_forest"],
            predictions_by_model["xgboost"],
            scaffolds,
            iterations=bootstrap_iterations,
            seed=seed + 7000,
        ) if {"xgboost", "random_forest"}.issubset(models) else None,
        "best_learned_vs_mean_baseline": _paired_cluster_bootstrap_difference(
            y,
            predictions_by_model["mean_baseline"],
            predictions_by_model[best],
            scaffolds,
            iterations=bootstrap_iterations,
            seed=seed + 8000,
        ),
    }
    best_oof = sorted((row for row in prediction_rows if row["model"] == best), key=lambda row: row["structure_id"])
    # Align by immutable structure identifier rather than relying on CSV order.
    prediction_by_structure = {row["structure_id"]: row for row in best_oof}
    residuals = np.asarray([
        abs(row["target_pactivity"] - prediction_by_structure[row["structure_id"]]["predicted_pactivity"])
        for row in rows
    ])
    final_half_width = _conformal_quantile(residuals)
    final_threshold = nearest_neighbor_threshold(packed)
    final_model = make_regressor(best, seed=seed)
    final_model.fit(X, y)
    bundle = {
        "artifact_version": 1,
        "model": final_model,
        "feature_config": config,
        "uncertainty": {
            "conformal_half_width_90": final_half_width,
            "coverage_target": 0.90,
            "source": "outer scaffold-held-out absolute residuals",
        },
        "applicability_domain": {
            "similarity_threshold": final_threshold,
            "reference_packed_fingerprints": packed,
            "method": "5th percentile training nearest-neighbor Morgan Tanimoto",
        },
        "training": {
            "model_name": best,
            "seed": seed,
            "target": dataset_metadata["assay_context"],
            "n_structures": len(rows),
            "n_scaffolds": len(unique_scaffolds),
            "target_definition": "median computed pActivity per standard InChIKey",
        },
    }
    model_path = output_dir / "best_model.joblib"
    joblib.dump(bundle, model_path, compress=3)
    _write_csv(output_dir / "fold_metrics.csv", fold_rows)
    _write_csv(output_dir / "oof_predictions.csv", prediction_rows)
    (output_dir / "feature_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    results = {
        "experiment": "PhytoTruth molecular bioactivity regression",
        "generated_at": datetime.now(UTC).isoformat(),
        "dataset": dataset_metadata,
        "validation": {
            "strategy": "5-fold Bemis-Murcko scaffold holdout with nested scaffold calibration split",
            "folds": folds,
            "unique_scaffolds": len(unique_scaffolds),
            "seed": seed,
            "bootstrap_iterations": bootstrap_iterations,
            "target_aggregation": "median per standardized InChIKey before splitting",
        },
        "models": models,
        "selection": {
            "primary_metric": "outer-fold MAE",
            "best_learned_model": best,
            "beats_mean_baseline_on_mae": best_mae < baseline_mae,
            "best_mae": best_mae,
            "mean_baseline_mae": baseline_mae,
        },
        "model_comparisons": comparisons,
        "serialized_model": {
            "file": model_path.name,
            "conformal_half_width_90": final_half_width,
            "applicability_similarity_threshold": final_threshold,
        },
        "evidence_reversal_status": "advantage_not_demonstrated; not used as a validated model feature",
    }
    results_path = output_dir / "metrics.json"
    results_path.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scikit_learn": sklearn.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "joblib": joblib.__version__,
        "rdkit": rdkit.__version__,
        "models_attempted": available_model_names(),
    }
    try:
        import xgboost

        environment["xgboost"] = xgboost.__version__
    except ImportError:
        environment["xgboost"] = None
    (output_dir / "environment.json").write_text(json.dumps(environment, indent=2, sort_keys=True), encoding="utf-8")
    _comparison_svg(output_dir / "model_comparison.svg", results)
    (output_dir / "model_card.md").write_text(_report(results), encoding="utf-8")
    command = (
        "python -m app.pipeline train-activity "
        f'--activities "{activities}" --natural-products "{natural_products}" '
        f'--output-dir "{output_dir}" --folds {folds} '
        f"--bootstrap-iterations {bootstrap_iterations} --seed {seed}"
    )
    (output_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    manifest = {
        "activities": str(activities.resolve()),
        "activities_sha256": _sha256(activities),
        "natural_products": str(natural_products.resolve()),
        "natural_products_sha256": _sha256(natural_products),
        "metrics_sha256": _sha256(results_path),
        "model_sha256": _sha256(model_path),
        "dataset_sha256": _sha256(output_dir / "training_dataset.csv"),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Train and validate PhytoTruth molecular bioactivity models")
    parser.add_argument("--activities", type=Path, required=True)
    parser.add_argument("--natural-products", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-key")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    results = run_training_experiment(
        args.activities,
        args.natural_products,
        args.output_dir,
        context_key=args.context_key,
        folds=args.folds,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "dataset": results["dataset"],
        "selection": results["selection"],
    }, indent=2))


if __name__ == "__main__":
    main()
