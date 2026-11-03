from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold

from .features import bemis_murcko_scaffold, feature_config, molecular_feature_vector
from .structures import standardize_smiles


SEED = 2026
POPCOUNT = np.asarray([int(value).bit_count() for value in range(256)], dtype=np.uint8)


@dataclass(frozen=True)
class FoldMetric:
    fold: int
    mae: float
    rmse: float
    r2: float
    train_size: int
    test_size: int


def regression_metrics(y_true: Sequence[float], prediction: Sequence[float]) -> dict[str, float]:
    truth = np.asarray(y_true, dtype=float)
    predicted = np.asarray(prediction, dtype=float)
    if truth.shape != predicted.shape or truth.ndim != 1 or not len(truth):
        raise ValueError("truth and prediction must be non-empty one-dimensional arrays of equal length")
    return {
        "mae": float(mean_absolute_error(truth, predicted)),
        "rmse": float(math.sqrt(mean_squared_error(truth, predicted))),
        "r2": float(r2_score(truth, predicted)) if len(truth) > 1 else float("nan"),
    }


def make_regressor(name: str, *, seed: int = SEED):
    if name == "mean_baseline":
        return DummyRegressor(strategy="mean")
    if name == "random_forest":
        return RandomForestRegressor(
            n_estimators=500,
            min_samples_leaf=2,
            max_features="sqrt",
            n_jobs=1,
            random_state=seed,
        )
    if name == "xgboost":
        try:
            from xgboost import XGBRegressor
        except ImportError as exc:  # pragma: no cover - depends on optional installation
            raise RuntimeError("XGBoost is not installed") from exc
        return XGBRegressor(
            n_estimators=500,
            max_depth=6,
            learning_rate=0.03,
            min_child_weight=3,
            subsample=0.8,
            colsample_bytree=0.7,
            reg_lambda=2.0,
            objective="reg:squarederror",
            tree_method="hist",
            n_jobs=1,
            random_state=seed,
        )
    raise ValueError(f"Unknown model: {name}")


def available_model_names() -> list[str]:
    names = ["mean_baseline", "random_forest"]
    try:
        import xgboost  # noqa: F401
    except ImportError:
        return names
    names.append("xgboost")
    return names


def scaffold_cross_validate(smiles: list[str], y: list[float], folds: int = 5) -> list[FoldMetric]:
    """Compatibility API: every standardized Murcko scaffold stays in one fold."""
    X = np.stack([molecular_feature_vector(item) for item in smiles])
    targets = np.asarray(y, dtype=float)
    groups = np.asarray([bemis_murcko_scaffold(item) for item in smiles])
    unique = np.unique(groups)
    if len(unique) < 2:
        raise ValueError("At least two unique scaffolds are required")
    splitter = GroupKFold(n_splits=min(folds, len(unique)))
    metrics: list[FoldMetric] = []
    for fold, (train, test) in enumerate(splitter.split(X, targets, groups), start=1):
        model = make_regressor("random_forest", seed=SEED + fold)
        model.fit(X[train], targets[train])
        values = regression_metrics(targets[test], model.predict(X[test]))
        metrics.append(FoldMetric(
            fold=fold,
            mae=values["mae"],
            rmse=values["rmse"],
            r2=values["r2"],
            train_size=len(train),
            test_size=len(test),
        ))
    return metrics


def pack_fingerprints(features: np.ndarray, n_bits: int = 2048) -> np.ndarray:
    matrix = np.asarray(features)
    if matrix.ndim != 2 or matrix.shape[1] < n_bits:
        raise ValueError("Feature matrix does not contain the configured fingerprint bits")
    bits = (matrix[:, :n_bits] > 0).astype(np.uint8)
    return np.packbits(bits, axis=1)


def max_tanimoto_similarity(query_packed: np.ndarray, reference_packed: np.ndarray, batch_size: int = 32) -> np.ndarray:
    query = np.asarray(query_packed, dtype=np.uint8)
    reference = np.asarray(reference_packed, dtype=np.uint8)
    if query.ndim != 2 or reference.ndim != 2 or query.shape[1] != reference.shape[1] or not len(reference):
        raise ValueError("Packed fingerprint matrices must be non-empty with matching widths")
    reference_counts = POPCOUNT[reference].sum(axis=1)
    maxima = np.empty(len(query), dtype=float)
    for start in range(0, len(query), batch_size):
        batch = query[start:start + batch_size]
        intersections = POPCOUNT[np.bitwise_and(batch[:, None, :], reference[None, :, :])].sum(axis=2)
        query_counts = POPCOUNT[batch].sum(axis=1)[:, None]
        unions = query_counts + reference_counts[None, :] - intersections
        similarities = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions > 0)
        maxima[start:start + len(batch)] = similarities.max(axis=1)
    return maxima


def nearest_neighbor_threshold(reference_packed: np.ndarray, quantile: float = 0.05, batch_size: int = 32) -> float:
    reference = np.asarray(reference_packed, dtype=np.uint8)
    if len(reference) < 2:
        raise ValueError("At least two reference molecules are required")
    counts = POPCOUNT[reference].sum(axis=1)
    nearest = np.empty(len(reference), dtype=float)
    for start in range(0, len(reference), batch_size):
        batch = reference[start:start + batch_size]
        intersections = POPCOUNT[np.bitwise_and(batch[:, None, :], reference[None, :, :])].sum(axis=2)
        unions = counts[start:start + len(batch), None] + counts[None, :] - intersections
        similarities = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions > 0)
        local_rows = np.arange(len(batch))
        similarities[local_rows, start + local_rows] = -1.0
        nearest[start:start + len(batch)] = similarities.max(axis=1)
    return float(np.quantile(nearest, quantile))


def load_model(path: Path | str) -> dict[str, Any]:
    import joblib

    bundle = joblib.load(path)
    required = {"model", "feature_config", "uncertainty", "applicability_domain", "training"}
    if not isinstance(bundle, dict) or not required.issubset(bundle):
        raise ValueError("Invalid PhytoTruth model bundle")
    return bundle


def predict_activity(bundle: dict[str, Any], smiles: str) -> dict[str, Any]:
    config = bundle["feature_config"]
    standardized = standardize_smiles(smiles)
    vector = molecular_feature_vector(
        standardized.canonical_smiles,
        radius=int(config["radius"]),
        n_bits=int(config["n_bits"]),
    )
    prediction = float(bundle["model"].predict(vector.reshape(1, -1))[0])
    packed = pack_fingerprints(vector.reshape(1, -1), int(config["n_bits"]))
    domain = bundle["applicability_domain"]
    similarity = float(max_tanimoto_similarity(packed, domain["reference_packed_fingerprints"])[0])
    threshold = float(domain["similarity_threshold"])
    base_half_width = float(bundle["uncertainty"]["conformal_half_width_90"])
    inflation = max(1.0, threshold / max(similarity, 1e-6))
    half_width = base_half_width * inflation
    return {
        "kind": "predicted",
        "predicted_pactivity": prediction,
        "prediction_interval_90": [prediction - half_width, prediction + half_width],
        "interval_method": "cross-validated absolute-residual conformal heuristic",
        "nearest_training_tanimoto": similarity,
        "applicability_domain_threshold": threshold,
        "inside_applicability_domain": similarity >= threshold,
        "canonical_smiles": standardized.canonical_smiles,
        "model_name": bundle["training"]["model_name"],
        "warning": None if similarity >= threshold else "Outside applicability domain; interval calibration is not guaranteed.",
    }
