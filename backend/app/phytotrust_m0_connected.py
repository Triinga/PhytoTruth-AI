from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from sklearn.metrics import mean_absolute_error, mean_squared_error

from .phytotrust_m0 import HierarchicalStudentTM0, M0Config, _student_scale


SEED = 2026
MODEL_VERSION = "phytotrust-m0-connected-v1"
EXPECTED_INPUT_NAME = "development_measurements.csv"
FORBIDDEN_INPUT_TOKENS = ("future_evaluation", "quarantined", "future_labels")
REQUIRED_COLUMNS = {
    "activity_id", "target_chembl_id", "molecule_id", "assay_id", "document_id",
    "document_year", "pchembl_value", "partition", "assay_context_key",
    "component_id", "is_verified_natural_product", "natural_product_sources",
}
METHODS = (
    "global_mean", "median_potency", "nonhierarchical_student_t",
    "phytotrust_m0_connected",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _joined(values: Iterable[Any]) -> str:
    return ";".join(sorted({str(value) for value in values if pd.notna(value)}))


def load_connected_development(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load only the public development artifact; fail closed for quarantined paths."""
    resolved = path.resolve()
    lowered = str(resolved).lower()
    if path.name.lower() != EXPECTED_INPUT_NAME or any(token in lowered for token in FORBIDDEN_INPUT_TOKENS):
        raise ValueError("M0-connected accepts only development_measurements.csv; evaluation artifacts are forbidden")
    frame = pd.read_csv(resolved)
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"Missing connected-M0 columns: {missing}")
    if frame.empty:
        raise ValueError("Connected development dataset is empty")
    if set(frame.partition.astype(str)) != {"development"}:
        raise AssertionError("Non-development rows entered connected M0")
    if (pd.to_numeric(frame.document_year, errors="coerce") >= 2020).any():
        raise AssertionError("Future evaluation years entered connected M0")
    if frame.activity_id.astype(str).duplicated().any():
        raise AssertionError("Duplicate activity identifiers in development input")
    if frame.target_chembl_id.nunique() != 1 or frame.target_chembl_id.iloc[0] != "CHEMBL222":
        raise AssertionError("Connected M0 is frozen to the retained CHEMBL222 context")
    if frame.assay_context_key.nunique() != 1:
        raise AssertionError("Incompatible assay contexts entered connected M0")
    if frame[["activity_id", "assay_id", "document_id", "molecule_id", "pchembl_value"]].isna().any().any():
        raise AssertionError("Required provenance or outcome is missing")
    frame = frame.copy()
    frame["document_year"] = frame.document_year.astype(int)
    frame["pchembl_value"] = frame.pchembl_value.astype(float)
    frame["is_verified_natural_product"] = frame.is_verified_natural_product.astype(bool)
    manifest = {
        "source_path": str(resolved),
        "source_sha256": _sha256(resolved),
        "raw_measurements": int(len(frame)),
        "compounds": int(frame.molecule_id.nunique()),
        "studies": int(frame.document_id.nunique()),
        "targets": sorted(frame.target_chembl_id.unique().tolist()),
        "contexts": sorted(frame.assay_context_key.unique().tolist()),
        "year_range": [int(frame.document_year.min()), int(frame.document_year.max())],
        "quarantined_files_opened": False,
    }
    return frame.sort_values(["document_year", "document_id", "molecule_id", "activity_id"]).reset_index(drop=True), manifest


def aggregate_compound_study(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse compatible within-publication records without combining publications."""
    rows: list[dict[str, Any]] = []
    keys = ["target_chembl_id", "assay_context_key", "component_id", "molecule_id", "document_id"]
    for key, group in frame.groupby(keys, sort=True, dropna=False):
        years = group.document_year.unique()
        if len(years) != 1:
            raise AssertionError(f"Publication has inconsistent years: {key[-1]}")
        rows.append({
            "target_chembl_id": key[0], "assay_context_key": key[1], "component_id": key[2],
            "molecule_id": key[3], "document_id": key[4], "document_year": int(years[0]),
            "study_median_pactivity": float(group.pchembl_value.median()),
            "within_study_min_pactivity": float(group.pchembl_value.min()),
            "within_study_max_pactivity": float(group.pchembl_value.max()),
            "raw_measurements": int(len(group)),
            "activity_ids": _joined(group.activity_id), "assay_ids": _joined(group.assay_id),
            "is_verified_natural_product": bool(group.is_verified_natural_product.any()),
            "natural_product_sources": _joined(
                source for value in group.natural_product_sources.fillna("")
                for source in str(value).split(";") if source
            ),
        })
    edges = pd.DataFrame(rows)
    if edges.duplicated(["molecule_id", "document_id"]).any():
        raise AssertionError("Compound-publication aggregation is not unique")
    return edges.sort_values(["document_year", "document_id", "molecule_id"]).reset_index(drop=True)


def construct_training_events(edges: pd.DataFrame) -> pd.DataFrame:
    """Create leave-publication-out residual events from genuine graph connections.

    Each offset excludes every record from the event's publication. This is used to
    estimate development variance, never to score the quarantined partition.
    """
    rows: list[dict[str, Any]] = []
    for molecule_id, compound in edges.groupby("molecule_id", sort=True):
        if compound.document_id.nunique() < 2:
            continue
        for heldout in compound.itertuples(index=False):
            prior = compound[compound.document_id != heldout.document_id]
            if prior.empty:
                continue
            rows.append({
                "target_chembl_id": heldout.target_chembl_id,
                "molecule_id": molecule_id,
                "heldout_study_key": heldout.document_id,
                "heldout_year": int(heldout.document_year),
                "prior_median_pactivity": float(prior.study_median_pactivity.median()),
                "heldout_median_pactivity": float(heldout.study_median_pactivity),
                "prior_study_count": int(prior.document_id.nunique()),
                "prior_document_ids": _joined(prior.document_id),
                "heldout_document_ids": str(heldout.document_id),
                "prior_assay_ids": _joined(";".join(prior.assay_ids).split(";")),
                "heldout_assay_ids": heldout.assay_ids,
                "prior_activity_ids": _joined(";".join(prior.activity_ids).split(";")),
                "heldout_activity_ids": heldout.activity_ids,
                "component_id": heldout.component_id,
                "is_verified_natural_product": bool(heldout.is_verified_natural_product),
                "history_mode": "leave_publication_out_development",
            })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(
        ["heldout_year", "heldout_study_key", "molecule_id"]
    ).reset_index(drop=True)


def construct_test_events(train_edges: pd.DataFrame, test_edges: pd.DataFrame, *, history_mode: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    train_by_compound = {key: value for key, value in train_edges.groupby("molecule_id", sort=False)}
    bridge_counts = test_edges[test_edges.molecule_id.isin(train_by_compound)].groupby("document_id").molecule_id.nunique()
    for heldout in test_edges.itertuples(index=False):
        prior = train_by_compound.get(heldout.molecule_id)
        if prior is None or prior.empty:
            continue
        if heldout.document_id in set(prior.document_id):
            raise AssertionError("Held-out publication entered historical evidence")
        prior_assays = set(";".join(prior.assay_ids).split(";"))
        prior_activities = set(";".join(prior.activity_ids).split(";"))
        if prior_assays & set(str(heldout.assay_ids).split(";")):
            raise AssertionError("Held-out assay entered historical evidence")
        if prior_activities & set(str(heldout.activity_ids).split(";")):
            raise AssertionError("Held-out activity entered historical evidence")
        count = int(prior.document_id.nunique())
        bridges = int(bridge_counts.get(heldout.document_id, 0))
        rows.append({
            "target_chembl_id": heldout.target_chembl_id,
            "molecule_id": heldout.molecule_id,
            "heldout_study_key": heldout.document_id,
            "heldout_year": int(heldout.document_year),
            "prior_median_pactivity": float(prior.study_median_pactivity.median()),
            "heldout_median_pactivity": float(heldout.study_median_pactivity),
            "prior_study_count": count,
            "prior_support_class": "multiple_prior_studies" if count >= 2 else "one_prior_study",
            "heldout_study_bridge_count": bridges,
            "study_bridge_class": "multi_bridge_study" if bridges >= 2 else "single_bridge_study",
            "prior_document_ids": _joined(prior.document_id),
            "heldout_document_ids": str(heldout.document_id),
            "prior_assay_ids": _joined(";".join(prior.assay_ids).split(";")),
            "heldout_assay_ids": heldout.assay_ids,
            "prior_activity_ids": _joined(";".join(prior.activity_ids).split(";")),
            "heldout_activity_ids": heldout.activity_ids,
            "component_id": heldout.component_id,
            "is_verified_natural_product": bool(heldout.is_verified_natural_product),
            "history_mode": history_mode,
        })
    return pd.DataFrame(rows)


def assert_fold_provenance(train_edges: pd.DataFrame, test_events: pd.DataFrame, *, temporal: bool) -> None:
    if test_events.empty:
        return
    heldout_docs = set(test_events.heldout_study_key)
    if set(train_edges.document_id) & heldout_docs:
        raise AssertionError("Publication leakage")
    train_assays = set(";".join(train_edges.assay_ids).split(";"))
    train_activities = set(";".join(train_edges.activity_ids).split(";"))
    test_assays = set(";".join(test_events.heldout_assay_ids).split(";"))
    test_activities = set(";".join(test_events.heldout_activity_ids).split(";"))
    if train_assays & test_assays:
        raise AssertionError("Assay leakage")
    if train_activities & test_activities:
        raise AssertionError("Activity leakage")
    for row in test_events.itertuples(index=False):
        if set(row.prior_document_ids.split(";")) & set(row.heldout_document_ids.split(";")):
            raise AssertionError("Compound history contains held-out publication")
        if set(row.prior_activity_ids.split(";")) & set(row.heldout_activity_ids.split(";")):
            raise AssertionError("Compound history contains held-out outcome")
    if temporal:
        if int(train_edges.document_year.max()) >= int(test_events.heldout_year.min()):
            raise AssertionError("Forward-temporal fold used contemporaneous or future evidence")


def _baseline_predictions(
    train_edges: pd.DataFrame, fit_events: pd.DataFrame, test: pd.DataFrame, config: M0Config,
) -> dict[str, dict[str, np.ndarray]]:
    outcome = test.heldout_median_pactivity.to_numpy(float)
    global_mean = float(train_edges.study_median_pactivity.mean())
    median_prediction = test.prior_median_pactivity.to_numpy(float)
    global_scale = _student_scale(
        train_edges.study_median_pactivity.to_numpy(float) - global_mean, config.degrees_of_freedom,
    )
    median_scale = _student_scale(
        fit_events.heldout_median_pactivity.to_numpy(float) - fit_events.prior_median_pactivity.to_numpy(float),
        config.degrees_of_freedom,
    )
    output: dict[str, dict[str, np.ndarray]] = {}
    for method, prediction, scale in (
        ("global_mean", np.full(len(test), global_mean), global_scale),
        ("median_potency", median_prediction, median_scale),
    ):
        values: dict[str, np.ndarray] = {
            "prediction": prediction,
            "nll": -student_t.logpdf((outcome - prediction) / scale, df=config.degrees_of_freedom) + math.log(scale),
            "residual_scale": np.full(len(test), scale),
            "between_study_sd": np.zeros(len(test)),
        }
        for level in (0.8, 0.9):
            critical = student_t.ppf((1.0 + level) / 2.0, df=config.degrees_of_freedom)
            suffix = str(int(level * 100))
            values[f"lower_{suffix}"] = prediction - critical * scale
            values[f"upper_{suffix}"] = prediction + critical * scale
        output[method] = values
    return output


def _fit_predict_fold(
    train_edges: pd.DataFrame, test_edges: pd.DataFrame, config: M0Config, *,
    scheme: str, fold: int, history_mode: str, temporal: bool,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    fit_events = construct_training_events(train_edges)
    test = construct_test_events(train_edges, test_edges, history_mode=history_mode)
    if len(fit_events) < 6 or fit_events.heldout_study_key.nunique() < 3 or test.empty:
        return pd.DataFrame(), []
    assert_fold_provenance(train_edges, test, temporal=temporal)
    nonhierarchical = HierarchicalStudentTM0(config, include_study_variance=False).fit(fit_events)
    m0 = HierarchicalStudentTM0(config, include_study_variance=True).fit(fit_events)
    predicted = {
        **_baseline_predictions(train_edges, fit_events, test, config),
        "nonhierarchical_student_t": nonhierarchical.predict_distribution(test, levels=(0.8, 0.9)),
        "phytotrust_m0_connected": m0.predict_distribution(test, levels=(0.8, 0.9)),
    }
    identity_columns = [
        "target_chembl_id", "molecule_id", "heldout_study_key", "heldout_year",
        "prior_median_pactivity", "heldout_median_pactivity", "prior_study_count",
        "prior_support_class", "heldout_study_bridge_count", "study_bridge_class",
        "prior_document_ids", "heldout_document_ids", "prior_assay_ids", "heldout_assay_ids",
        "prior_activity_ids", "heldout_activity_ids", "component_id",
        "is_verified_natural_product", "history_mode",
    ]
    identity = test[identity_columns].copy().reset_index(drop=True)
    identity.insert(0, "prediction_id", [f"{scheme}:{fold}:{i}" for i in range(len(identity))])
    rows = []
    for method, values in predicted.items():
        result = identity.copy()
        result.insert(0, "method", method)
        result.insert(0, "fold", fold)
        result.insert(0, "validation_scheme", scheme)
        for name, value in values.items():
            result[name] = value
        rows.append(result)
    parameter_rows = []
    for method, model in (("nonhierarchical_student_t", nonhierarchical), ("phytotrust_m0_connected", m0)):
        parameters = model.parameters()
        parameter_rows.append({
            "validation_scheme": scheme, "fold": fold, "method": method,
            "training_edges": int(len(train_edges)), "training_fit_events": int(len(fit_events)),
            "training_connected_compounds": int(fit_events.molecule_id.nunique()),
            "training_connected_studies": int(fit_events.heldout_study_key.nunique()),
            "test_events": int(len(test)), "test_studies": int(test.heldout_study_key.nunique()),
            "alpha": float(parameters["target_effects"]["CHEMBL222"]),
            "residual_scale": float(parameters["residual_scale"]),
            "between_study_sd": float(parameters["target_between_study_sd"]["CHEMBL222"]),
            "optimizer_success": bool(parameters["optimization"]["success"]),
        })
    return pd.concat(rows, ignore_index=True), parameter_rows


def cross_validate(edges: pd.DataFrame, config: M0Config) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    predictions: list[pd.DataFrame] = []
    parameters: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    repeated = edges.groupby("molecule_id").document_id.transform("nunique") >= 2
    eligible_studies = sorted(edges.loc[repeated, "document_id"].unique())
    for fold, study in enumerate(eligible_studies):
        train = edges[edges.document_id != study].copy()
        test = edges[edges.document_id == study].copy()
        result, params = _fit_predict_fold(
            train, test, config, scheme="study_grouped", fold=fold,
            history_mode="other_development_studies_not_time_restricted", temporal=False,
        )
        if result.empty:
            skipped.append({"scheme": "study_grouped", "fold": fold, "group": study, "reason": "insufficient_connected_training_data"})
        else:
            predictions.append(result)
            parameters.extend(params)
    eligible_years = sorted(edges.document_year.unique())
    temporal_fold = 0
    for year in eligible_years:
        train = edges[edges.document_year < year].copy()
        test = edges[edges.document_year == year].copy()
        if train.empty:
            continue
        result, params = _fit_predict_fold(
            train, test, config, scheme="forward_temporal", fold=temporal_fold,
            history_mode="strictly_earlier_publications", temporal=True,
        )
        if result.empty:
            skipped.append({"scheme": "forward_temporal", "fold": temporal_fold, "group": int(year), "reason": "no_eligible_transfer_events_or_connected_training_data"})
        else:
            predictions.append(result)
            parameters.extend(params)
        temporal_fold += 1
    if not predictions:
        raise RuntimeError("No connected-M0 validation folds were estimable")
    return pd.concat(predictions, ignore_index=True), pd.DataFrame(parameters), {"skipped_folds": skipped}


def _metric_row(group: pd.DataFrame) -> dict[str, Any]:
    observed = group.heldout_median_pactivity.to_numpy(float)
    predicted = group.prediction.to_numpy(float)
    result: dict[str, Any] = {
        "n": int(len(group)), "compounds": int(group.molecule_id.nunique()),
        "studies": int(group.heldout_study_key.nunique()),
        "mae": float(mean_absolute_error(observed, predicted)),
        "rmse": float(math.sqrt(mean_squared_error(observed, predicted))),
        "nll": float(group.nll.mean()),
    }
    for level in (80, 90):
        covered = (observed >= group[f"lower_{level}"]) & (observed <= group[f"upper_{level}"])
        result[f"coverage_{level}"] = float(covered.mean())
        result[f"mean_width_{level}"] = float((group[f"upper_{level}"] - group[f"lower_{level}"]).mean())
    return result


def evaluate_predictions(predictions: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}
    scopes = {
        "all": lambda value: value,
        "natural_products": lambda value: value[value.is_verified_natural_product],
        "one_prior_study": lambda value: value[value.prior_support_class == "one_prior_study"],
        "multiple_prior_studies": lambda value: value[value.prior_support_class == "multiple_prior_studies"],
        "single_bridge_study": lambda value: value[value.study_bridge_class == "single_bridge_study"],
        "multi_bridge_study": lambda value: value[value.study_bridge_class == "multi_bridge_study"],
    }
    for (scheme, method), group in predictions.groupby(["validation_scheme", "method"], sort=True):
        for scope, selector in scopes.items():
            scoped = selector(group)
            values = _metric_row(scoped) if len(scoped) else {}
            rows.append({"validation_scheme": scheme, "scope": scope, "method": method, **values})
            nested.setdefault(scheme, {}).setdefault(scope, {})[method] = values
    return pd.DataFrame(rows), nested


def study_cluster_bootstrap(
    predictions: pd.DataFrame, *, iterations: int, seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}
    for scheme in sorted(predictions.validation_scheme.unique()):
        scheme_data = predictions[predictions.validation_scheme == scheme]
        for scope, scoped in (
            ("all", scheme_data),
            ("natural_products", scheme_data[scheme_data.is_verified_natural_product]),
        ):
            studies = sorted(scoped.heldout_study_key.unique())
            if len(studies) < 2:
                continue
            nested.setdefault(scheme, {}).setdefault(scope, {})
            by_method = {method: scoped[scoped.method == method].set_index("prediction_id") for method in METHODS}
            common = set.intersection(*(set(value.index) for value in by_method.values()))
            if not common:
                continue
            order = sorted(common)
            base = by_method["phytotrust_m0_connected"].loc[order]
            clusters = {study: base.index[base.heldout_study_key == study].tolist() for study in studies}
            observed = base.heldout_median_pactivity.to_numpy(float)
            for comparator in METHODS[:-1]:
                comparison = by_method[comparator].loc[order]
                point = {
                    "mae": float(np.mean(np.abs(observed - base.prediction)) - np.mean(np.abs(observed - comparison.prediction))),
                    "rmse": float(np.sqrt(np.mean((observed - base.prediction) ** 2)) - np.sqrt(np.mean((observed - comparison.prediction) ** 2))),
                    "nll": float(base.nll.mean() - comparison.nll.mean()),
                }
                samples = {metric: [] for metric in point}
                for _ in range(iterations):
                    chosen = rng.choice(studies, len(studies), replace=True)
                    ids = [item for study in chosen for item in clusters[study]]
                    m0_sample = base.loc[ids]
                    other_sample = comparison.loc[ids]
                    truth = m0_sample.heldout_median_pactivity.to_numpy(float)
                    samples["mae"].append(float(np.mean(np.abs(truth - m0_sample.prediction)) - np.mean(np.abs(truth - other_sample.prediction))))
                    samples["rmse"].append(float(np.sqrt(np.mean((truth - m0_sample.prediction) ** 2)) - np.sqrt(np.mean((truth - other_sample.prediction) ** 2))))
                    samples["nll"].append(float(m0_sample.nll.mean() - other_sample.nll.mean()))
                ci = {metric: [float(value) for value in np.quantile(sample, [0.025, 0.975])] for metric, sample in samples.items()}
                significant = {metric: bool(interval[1] < 0 or interval[0] > 0) for metric, interval in ci.items()}
                nested[scheme][scope][comparator] = {"point": point, "ci95": ci, "ci_excludes_zero": significant}
                records.append({
                    "validation_scheme": scheme, "scope": scope,
                    "comparison": f"phytotrust_m0_connected-minus-{comparator}",
                    **{f"delta_{key}": value for key, value in point.items()},
                    **{f"delta_{key}_ci95": json.dumps(value) for key, value in ci.items()},
                    **{f"delta_{key}_ci_excludes_zero": value for key, value in significant.items()},
                })
    return pd.DataFrame(records), nested


def parameter_bootstrap(
    fit_events: pd.DataFrame, fitted: HierarchicalStudentTM0, *, iterations: int, seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    params = fitted.parameters()
    alpha = float(params["target_effects"]["CHEMBL222"])
    sigma = float(params["residual_scale"])
    tau = float(params["target_between_study_sd"]["CHEMBL222"])
    studies = sorted(fit_events.heldout_study_key.unique())
    records = []
    for iteration in range(iterations):
        shifts = dict(zip(studies, rng.normal(0.0, tau, len(studies))))
        simulated = fit_events.copy()
        simulated["heldout_median_pactivity"] = (
            simulated.prior_median_pactivity.to_numpy(float) + alpha
            + simulated.heldout_study_key.map(shifts).to_numpy(float)
            + rng.standard_t(fitted.config.degrees_of_freedom, len(simulated)) * sigma
        )
        try:
            model = HierarchicalStudentTM0(fitted.config, include_study_variance=True).fit(simulated)
            estimate = model.parameters()
            records.append({
                "iteration": iteration,
                "alpha": float(estimate["target_effects"]["CHEMBL222"]),
                "residual_scale": float(estimate["residual_scale"]),
                "between_study_sd": float(estimate["target_between_study_sd"]["CHEMBL222"]),
                "optimizer_success": bool(estimate["optimization"]["success"]),
            })
        except (RuntimeError, ValueError, OverflowError):
            records.append({"iteration": iteration, "optimizer_success": False})
    table = pd.DataFrame(records)
    successful = table[table.optimizer_success & table.between_study_sd.notna()]
    if successful.empty:
        summary = {"successful": 0, "requested": iterations, "identifiable": False}
    else:
        summary = {
            "requested": iterations, "successful": int(len(successful)),
            "between_study_sd_point": tau,
            "between_study_sd_ci95_model_based": [float(value) for value in successful.between_study_sd.quantile([0.025, 0.975])],
            "residual_scale_ci95_model_based": [float(value) for value in successful.residual_scale.quantile([0.025, 0.975])],
            "boundary_fraction_tau_below_0_01": float((successful.between_study_sd < 0.01).mean()),
            "note": "Parametric bootstrap is conditional on the fitted M0 and does not correct graph-selection or publication bias.",
        }
    return table, summary


def identifiability_diagnostics(
    edges: pd.DataFrame, fold_parameters: pd.DataFrame, final_model: HierarchicalStudentTM0,
    parameter_uncertainty: dict[str, Any],
) -> dict[str, Any]:
    counts = edges.groupby("molecule_id").document_id.nunique()
    repeated_ids = counts[counts >= 2].index
    connected = edges[edges.molecule_id.isin(repeated_ids)]
    study_anchor_counts = connected.groupby("document_id").molecule_id.nunique()
    tau_folds = fold_parameters.loc[
        fold_parameters.method == "phytotrust_m0_connected", "between_study_sd"
    ].dropna()
    final = final_model.parameters()
    structural_gate = bool(
        len(repeated_ids) >= 20
        and int((study_anchor_counts >= 2).sum()) >= 10
        and int((study_anchor_counts >= 3).sum()) >= 5
    )
    bootstrap_ci = parameter_uncertainty.get("between_study_sd_ci95_model_based", [0.0, float("inf")])
    stable = bool(len(tau_folds) and tau_folds.min() >= 0.01 and tau_folds.max() / max(tau_folds.min(), 1e-9) <= 4.0)
    uncertainty_away_from_boundary = bool(bootstrap_ci[0] >= 0.01)
    identifiable = structural_gate and stable and uncertainty_away_from_boundary
    return {
        "criterion": (
            ">=20 repeated compounds, >=10 studies with >=2 shared-compound anchors, "
            ">=5 with >=3 anchors, fold tau range <=4-fold without boundary collapse, and model-based 95% interval above 0.01"
        ),
        "repeated_compounds": int(len(repeated_ids)),
        "connected_edges": int(len(connected)),
        "connected_studies": int(connected.document_id.nunique()),
        "studies_with_at_least_2_shared_compounds": int((study_anchor_counts >= 2).sum()),
        "studies_with_at_least_3_shared_compounds": int((study_anchor_counts >= 3).sum()),
        "compounds_in_at_least_3_studies": int((counts >= 3).sum()),
        "structural_gate_passed": structural_gate,
        "final_between_study_sd": float(final["target_between_study_sd"]["CHEMBL222"]),
        "cross_fold_tau_range": [float(tau_folds.min()), float(tau_folds.max())] if len(tau_folds) else None,
        "fold_stability_passed": stable,
        "model_based_uncertainty_away_from_boundary": uncertainty_away_from_boundary,
        "between_study_variance_identifiable": identifiable,
        "conclusion": (
            "Between-study variance met all prespecified development-only identifiability checks. External validation is still required."
            if identifiable else
            "Between-study variance is not reliably identifiable: the connected graph and/or fold/bootstrap stability criteria failed. The estimate is a regularized feasibility quantity, not a definitive variance component."
        ),
    }


def _render_report(results: dict[str, Any], metrics: pd.DataFrame) -> str:
    dataset = results["dataset"]
    ident = results["identifiability"]
    lines = [
        "# PhytoTRUST-M0 on the study-connected CHEMBL222 dataset", "",
        "## Scope and leakage contract", "",
        f"This experiment used only `{EXPECTED_INPUT_NAME}`: {dataset['raw_measurements']} real development measurements, "
        f"{dataset['compounds']} compounds, and {dataset['studies']} ChEMBL publications. No future-evaluation history or quarantined label file was opened.", "",
        "Study-grouped validation leaves out entire publications and excludes their activity and assay identifiers. It is a retrospective domain-holdout analysis and is not described as prospective because other development publications can be later in time. Forward-temporal validation uses only strictly earlier publication years. Same-compound historical measurements are intentionally allowed because prediction for a previously measured compound is the declared estimand; the held-out publication's records are never used in its features or fitting data.", "",
        "M0 uses the prespecified Student-t residual model with a normally distributed publication random intercept. The random effect for a new publication is integrated out. The non-hierarchical ablation has the same location and Student-t likelihood but fixes between-study variance to zero.", "",
        "## Validation metrics", "",
        "| Validation | Scope | Method | n | studies | MAE | RMSE | NLL | 80% coverage / width | 90% coverage / width |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    shown = metrics[metrics.scope.isin(["all", "natural_products", "one_prior_study", "multiple_prior_studies"])]
    for row in shown.itertuples():
        if pd.isna(getattr(row, "n", np.nan)):
            continue
        lines.append(
            f"| {row.validation_scheme} | {row.scope} | {row.method} | {int(row.n)} | {int(row.studies)} | "
            f"{row.mae:.3f} | {row.rmse:.3f} | {row.nll:.3f} | "
            f"{row.coverage_80:.3f} / {row.mean_width_80:.3f} | {row.coverage_90:.3f} / {row.mean_width_90:.3f} |"
        )
    sg = results["metrics"].get("study_grouped", {}).get("all", {})
    ft = results["metrics"].get("forward_temporal", {}).get("all", {})
    m0_sg = sg.get("phytotrust_m0_connected", {})
    m0_ft = ft.get("phytotrust_m0_connected", {})
    median_sg = sg.get("median_potency", {})
    median_ft = ft.get("median_potency", {})
    nonhier_sg = sg.get("nonhierarchical_student_t", {})
    comparison = results["study_cluster_bootstrap_differences"]
    sg_ci = comparison.get("study_grouped", {}).get("all", {}).get("median_potency", {}).get("ci95", {})
    sensitivity = results["dominant_study_sensitivity"]
    lines += [
        "", "## Primary comparison", "",
        f"Study-grouped M0: MAE {m0_sg.get('mae', float('nan')):.3f}, RMSE {m0_sg.get('rmse', float('nan')):.3f}, NLL {m0_sg.get('nll', float('nan')):.3f}. Median potency: MAE {median_sg.get('mae', float('nan')):.3f}, RMSE {median_sg.get('rmse', float('nan')):.3f}, NLL {median_sg.get('nll', float('nan')):.3f}. Matched no-study-variance ablation NLL: {nonhier_sg.get('nll', float('nan')):.3f}.", "",
        f"Forward-temporal M0: MAE {m0_ft.get('mae', float('nan')):.3f}, RMSE {m0_ft.get('rmse', float('nan')):.3f}, NLL {m0_ft.get('nll', float('nan')):.3f}. Median potency: MAE {median_ft.get('mae', float('nan')):.3f}, RMSE {median_ft.get('rmse', float('nan')):.3f}, NLL {median_ft.get('nll', float('nan')):.3f}.", "",
        f"The study-cluster bootstrap 95% interval for study-grouped M0-minus-median NLL is {sg_ci.get('nll', 'not estimable')}. Negative values favor M0. These intervals use publication as the resampling unit and are not independent confirmation.", "",
        "## Ablations and dominant-publication sensitivity", "",
        "The matched non-hierarchical Student-t ablation removes only the publication random-effect variance. Median potency removes the learned global correction and uses a single empirical residual scale. Global mean removes compound history as well. M0 was not significantly better than median potency or the matched no-study-variance ablation in the publication-cluster bootstrap.", "",
        f"The largest connected publication was {sensitivity['removed_study']} ({sensitivity['removed_transfer_events']} held-out transfer events). After excluding those held-out predictions, study-grouped M0 NLL was {sensitivity['metrics_without_study']['phytotrust_m0_connected']['nll']:.3f}, versus {sensitivity['metrics_without_study']['nonhierarchical_student_t']['nll']:.3f} for the no-study-variance ablation. Refitting the final variance model without that publication changed between-study SD from {ident['final_between_study_sd']:.4f} to {sensitivity['refit_between_study_sd']:.4f} pChEMBL.", "",
        "## Connectivity analysis", "",
        f"The variance-estimation subgraph contains {ident['connected_edges']} edges from {ident['repeated_compounds']} repeated compounds and {ident['connected_studies']} publications. Only {ident['studies_with_at_least_2_shared_compounds']} publications have at least two shared-compound anchors, {ident['studies_with_at_least_3_shared_compounds']} have at least three, and {ident['compounds_in_at_least_3_studies']} compound appears in at least three studies.", "",
        "Metrics stratified by one versus multiple prior studies and by single- versus multi-bridge publications are included in `metrics.csv`. The multiple-prior stratum is too small for a stable claim; differences are descriptive rather than causal evidence that connectivity improves prediction.", "",
        "## Parameter uncertainty and identifiability", "",
        ident["conclusion"], "",
        f"Final between-study SD: {ident['final_between_study_sd']:.4f} pChEMBL; cross-fold range: {ident['cross_fold_tau_range']}. Model-based parametric-bootstrap 95% interval: {results['parameter_uncertainty'].get('between_study_sd_ci95_model_based', 'not estimable')}. The bootstrap is conditional on M0 and cannot account for publication selection, dependence induced by leave-publication-out medians, or database curation uncertainty.", "",
        "## Natural products", "",
        f"Only {dataset['verified_natural_product_compounds']} verified natural-product compounds occur in the full development input, and {dataset['verified_natural_product_transfer_events']} appear in at least one validation event. NP status is used only for subgroup reporting, never as a predictor or fitting weight. Subgroup estimates are exploratory.", "",
        "## Limitations", "",
        "- M0 predicts a comparable experimental measurement, not true potency, efficacy, safety, or clinical benefit.",
        "- Publication is a proxy for study/laboratory and may merge or split actual experimental domains.",
        "- Leave-publication-out residual events share compound measurements and are statistically dependent; the study-cluster bootstrap does not remove all crossed dependence.",
        "- Study-grouped validation is not temporally prospective; only the forward-temporal analysis enforces the time filtration.",
        "- Sparse repeated compounds, a single target/context, and very few natural products limit generalization.",
        "- A structural graph gate is necessary but not sufficient for variance identifiability.",
        "", "## Scientific verdict", "",
        results["scientific_verdict"],
    ]
    return "\n".join(lines) + "\n"


def _plot(metrics: pd.DataFrame, path: Path) -> None:
    rows = metrics[(metrics.scope == "all") & metrics.method.isin(METHODS)]
    colors = {"global_mean": "#64748b", "median_potency": "#0ea5e9", "nonhierarchical_student_t": "#a78bfa", "phytotrust_m0_connected": "#22c55e"}
    labels = {"global_mean": "Global mean", "median_potency": "Median potency", "nonhierarchical_student_t": "No study variance", "phytotrust_m0_connected": "PhytoTRUST-M0"}
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="640" viewBox="0 0 1120 640">',
        '<rect width="1120" height="640" fill="#07111f"/>',
        '<text x="55" y="42" fill="#f8fafc" font-family="sans-serif" font-size="22">CHEMBL222 connected-data M0 validation</text>',
        '<text x="55" y="66" fill="#94a3b8" font-family="sans-serif" font-size="13">Development data only; lower is better. Error bars are reported in CSV/JSON.</text>',
    ]
    for panel, scheme in enumerate(("study_grouped", "forward_temporal")):
        subset = rows[rows.validation_scheme == scheme]
        x0 = 55 + 535 * panel
        pieces.append(f'<text x="{x0}" y="105" fill="#e2e8f0" font-family="sans-serif" font-size="17">{scheme.replace("_", " ").title()}</text>')
        for metric_index, metric in enumerate(("nll", "mae")):
            y0 = 140 + 225 * metric_index
            maximum = max(float(subset[metric].max()), 0.001)
            pieces.append(f'<text x="{x0}" y="{y0}" fill="#94a3b8" font-family="sans-serif" font-size="13">{metric.upper()}</text>')
            for index, method in enumerate(METHODS):
                match = subset[subset.method == method]
                if match.empty:
                    continue
                value = float(match.iloc[0][metric])
                y = y0 + 24 + index * 42
                width = 290 * value / maximum
                pieces.append(f'<text x="{x0}" y="{y+15}" fill="#cbd5e1" font-family="sans-serif" font-size="12">{labels[method]}</text>')
                pieces.append(f'<rect x="{x0+145}" y="{y}" width="{width:.1f}" height="22" rx="4" fill="{colors[method]}"/>')
                pieces.append(f'<text x="{x0+151+width:.1f}" y="{y+15}" fill="#f8fafc" font-family="sans-serif" font-size="12">{value:.3f}</text>')
    pieces.append("</svg>")
    path.write_text("\n".join(pieces), encoding="utf-8")


def run_experiment(
    dataset: Path, output_dir: Path, *, bootstrap_iterations: int = 1000,
    parameter_bootstrap_iterations: int = 200, seed: int = SEED,
    config: M0Config | None = None,
) -> dict[str, Any]:
    config = config or M0Config(quadrature_nodes=15, optimizer_maxiter=700)
    frame, manifest = load_connected_development(dataset)
    edges = aggregate_compound_study(frame)
    predictions, fold_parameters, split_details = cross_validate(edges, config)
    metrics_table, metrics = evaluate_predictions(predictions)
    bootstrap_table, bootstrap = study_cluster_bootstrap(
        predictions, iterations=bootstrap_iterations, seed=seed,
    )
    final_events = construct_training_events(edges)
    final_model = HierarchicalStudentTM0(config, include_study_variance=True).fit(final_events)
    parameter_table, parameter_uncertainty = parameter_bootstrap(
        final_events, final_model, iterations=parameter_bootstrap_iterations, seed=seed + 1,
    )
    identifiability = identifiability_diagnostics(edges, fold_parameters, final_model, parameter_uncertainty)
    repeated_ids = set(edges.groupby("molecule_id").document_id.nunique().loc[lambda x: x >= 2].index)
    manifest.update({
        "compound_study_edges": int(len(edges)),
        "connected_training_events": int(len(final_events)),
        "repeated_compounds": int(len(repeated_ids)),
        "connected_studies": int(final_events.heldout_study_key.nunique()),
        "verified_natural_product_compounds": int(edges.loc[edges.is_verified_natural_product, "molecule_id"].nunique()),
        "verified_natural_product_transfer_events": int(final_events.is_verified_natural_product.sum()),
    })
    largest_connected_study = str(final_events.heldout_study_key.value_counts().index[0])
    sensitivity_predictions = predictions[
        ~((predictions.validation_scheme == "study_grouped")
          & (predictions.heldout_study_key == largest_connected_study))
    ]
    sensitivity_table, sensitivity_nested = evaluate_predictions(sensitivity_predictions)
    sensitivity_all = sensitivity_nested["study_grouped"]["all"]
    sensitivity_edges = edges[edges.document_id != largest_connected_study]
    sensitivity_events = construct_training_events(sensitivity_edges)
    sensitivity_model = HierarchicalStudentTM0(config, include_study_variance=True).fit(sensitivity_events)
    sensitivity_tau = float(sensitivity_model.parameters()["target_between_study_sd"]["CHEMBL222"])
    sg = metrics.get("study_grouped", {}).get("all", {})
    ft = metrics.get("forward_temporal", {}).get("all", {})
    sg_delta = bootstrap.get("study_grouped", {}).get("all", {}).get("median_potency", {})
    supported = bool(
        identifiability["between_study_variance_identifiable"]
        and sg_delta.get("ci95", {}).get("nll", [0, 0])[1] < 0
        and ft.get("phytotrust_m0_connected", {}).get("nll", math.inf)
        < ft.get("median_potency", {}).get("nll", -math.inf)
    )
    verdict = (
        "Development-only evidence supports M0 as a variance-aware foundation, subject to a new untouched external evaluation. It does not establish methodological novelty."
        if supported else
        "Development-only evidence does not establish M0 as a superior or statistically identified foundation. The connected graph permits a more relevant feasibility analysis, but sparse anchors and/or inconsistent held-out performance require more connected independent studies before M1 or a journal-level claim."
    )
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "model_version": MODEL_VERSION,
        "seed": seed, "configuration": asdict(config), "dataset": manifest,
        "validation": {
            "study_grouped": "leave-one-publication-out; retrospective, not time restricted",
            "forward_temporal": "train publication year strictly earlier than test year",
            "heldout_publication_assay_activity_excluded": True,
            "same_compound_history_allowed_by_estimand": True,
            "quarantined_artifacts_loaded": False,
            **split_details,
        },
        "metrics": metrics,
        "study_cluster_bootstrap_differences": bootstrap,
        "parameter_uncertainty": parameter_uncertainty,
        "identifiability": identifiability,
        "dominant_study_sensitivity": {
            "removed_study": largest_connected_study,
            "removed_transfer_events": int((
                (predictions.validation_scheme == "study_grouped")
                & (predictions.method == "phytotrust_m0_connected")
                & (predictions.heldout_study_key == largest_connected_study)
            ).sum()),
            "refit_between_study_sd": sensitivity_tau,
            "metrics_without_study": sensitivity_all,
        },
        "final_parameters": final_model.parameters(),
        "supported_as_foundation": supported,
        "scientific_verdict": verdict,
        "scientific_status": "development_only_connected_m0_feasibility",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "predictions.csv", index=False)
    metrics_table.to_csv(output_dir / "metrics.csv", index=False)
    fold_parameters.to_csv(output_dir / "fold_parameters.csv", index=False)
    bootstrap_table.to_csv(output_dir / "bootstrap_differences.csv", index=False)
    parameter_table.to_csv(output_dir / "parameter_bootstrap.csv", index=False)
    sensitivity_table.to_csv(output_dir / "dominant_study_sensitivity.csv", index=False)
    final_events.to_csv(output_dir / "connected_variance_events.csv", index=False)
    edges.to_csv(output_dir / "compound_study_aggregates.csv", index=False)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir / "model_config.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "configuration": asdict(config),
        "equation": "Y_cs = median(Y_c,other allowed studies) + alpha + u_s + StudentT(df, 0, sigma)",
        "study_effect": "u_s ~ Normal(0, tau^2); u_new integrated out for blind transfer",
        "quarantine_contract": "Only development_measurements.csv may be loaded",
    }, indent=2), encoding="utf-8")
    (output_dir / "model_parameters.json").write_text(json.dumps(final_model.parameters(), indent=2), encoding="utf-8")
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_render_report(results, metrics_table), encoding="utf-8")
    _plot(metrics_table, output_dir / "comparison.svg")
    command = (
        "python -m app.phytotrust_m0_connected "
        "--dataset ../experiments/phytotrust_connected_dataset/development_measurements.csv "
        "--output-dir ../experiments/phytotrust_m0_connected "
        f"--bootstrap-iterations {bootstrap_iterations} "
        f"--parameter-bootstrap-iterations {parameter_bootstrap_iterations} --seed {seed}\n"
    )
    (output_dir / "command.txt").write_text(command, encoding="utf-8")
    artifact_names = sorted(path.name for path in output_dir.iterdir() if path.is_file() and path.name != "manifest.json")
    artifact_hashes = {name: _sha256(output_dir / name) for name in artifact_names}
    (output_dir / "manifest.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "source_dataset": str(dataset.resolve()),
        "source_sha256": manifest["source_sha256"], "files": artifact_hashes,
        "quarantined_artifacts_loaded": False,
    }, indent=2), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="PhytoTRUST-M0 on connected CHEMBL222 development data")
    parser.add_argument("--dataset", type=Path, default=Path("../experiments/phytotrust_connected_dataset/development_measurements.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("../experiments/phytotrust_m0_connected"))
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--parameter-bootstrap-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    results = run_experiment(
        args.dataset, args.output_dir, bootstrap_iterations=args.bootstrap_iterations,
        parameter_bootstrap_iterations=args.parameter_bootstrap_iterations, seed=args.seed,
    )
    print(json.dumps({
        "dataset": results["dataset"], "metrics": results["metrics"],
        "identifiability": results["identifiability"],
        "supported_as_foundation": results["supported_as_foundation"],
    }, indent=2))


if __name__ == "__main__":
    main()
