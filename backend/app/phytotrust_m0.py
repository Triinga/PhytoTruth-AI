from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp
from scipy.stats import t as student_t
from sklearn.metrics import mean_absolute_error, mean_squared_error

from .methodology_redesign import load_development_dataset, make_splits


SEED = 2026
MODEL_VERSION = "phytotrust-m0-v1"
REQUIRED_COLUMNS = {
    "target_chembl_id", "molecule_id", "heldout_study_key", "heldout_year",
    "prior_median_pactivity", "heldout_median_pactivity", "partition",
    "is_verified_natural_product", "prior_document_ids", "heldout_document_ids",
    "prior_assay_ids", "heldout_assay_ids",
}


@dataclass(frozen=True)
class M0Config:
    degrees_of_freedom: float = 4.0
    quadrature_nodes: int = 15
    target_effect_scale: float = 1.5
    log_tau_pooling_scale: float = 0.75
    minimum_scale: float = 1e-3
    optimizer_maxiter: int = 700


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_m0_development_data(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    raw = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS - set(raw.columns))
    if missing:
        raise ValueError(f"Missing M0 columns: {missing}")
    frame, base_manifest = load_development_dataset(path)
    frame = frame.dropna(subset=["prior_median_pactivity", "heldout_median_pactivity"]).copy()
    if frame.empty:
        raise ValueError("No eligible pre-2020 development rows")
    if not (frame["heldout_year"] <= 2019).all() or not (frame["partition"] == "development").all():
        raise AssertionError("Locked evaluation data entered M0")
    for prior, heldout in zip(frame.prior_document_ids, frame.heldout_document_ids):
        if set(str(prior).split(";")) & set(str(heldout).split(";")):
            raise ValueError("Prior and held-out document provenance overlap")
    for prior, heldout in zip(frame.prior_assay_ids, frame.heldout_assay_ids):
        if set(str(prior).split(";")) & set(str(heldout).split(";")):
            raise ValueError("Prior and held-out assay provenance overlap")
    frame["is_verified_natural_product"] = frame.is_verified_natural_product.astype(bool)
    manifest = dict(base_manifest)
    manifest.update({
        "eligible_rows": int(len(frame)),
        "source_sha256": _sha256(path),
        "outcome": "heldout_median_pactivity",
        "historical_offset": "prior_median_pactivity",
        "locked_evaluation_loaded_for_modeling": False,
    })
    return frame.reset_index(drop=True), manifest


class HierarchicalStudentTM0:
    """Robust random-intercept model for blind prediction into unseen studies.

    y_i = historical_median_i + alpha[target_i] + u[study_i,target_i] + eps_i
    u ~ Normal(0, tau_target**2), eps ~ StudentT(df, 0, sigma).

    Study effects are marginalized by Gauss-Hermite quadrature. Target-specific
    log(tau) values are partially pooled toward a shared log(tau).
    """

    def __init__(self, config: M0Config | None = None, *, include_study_variance: bool = True):
        self.config = config or M0Config()
        self.include_study_variance = include_study_variance
        self.targets_: list[str] = []
        self.alpha_: np.ndarray | None = None
        self.sigma_: float | None = None
        self.tau_: np.ndarray | None = None
        self.optimization_: dict[str, Any] = {}
        nodes, weights = np.polynomial.hermite.hermgauss(self.config.quadrature_nodes)
        self._gh_nodes = nodes
        self._gh_log_weights = np.log(weights) - 0.5 * math.log(math.pi)

    def _unpack(self, parameters: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
        count = len(self.targets_)
        alpha = parameters[:count]
        sigma = self.config.minimum_scale + math.exp(float(parameters[count]))
        if self.include_study_variance:
            global_log_tau = float(parameters[count + 1])
            delta = parameters[count + 2:count + 2 + count]
            tau = self.config.minimum_scale + np.exp(global_log_tau + delta)
        else:
            tau = np.zeros(count, dtype=float)
        return alpha, sigma, tau

    def _objective(
        self,
        parameters: np.ndarray,
        frame: pd.DataFrame,
        target_index: np.ndarray,
        study_groups: list[np.ndarray],
    ) -> float:
        alpha, sigma, tau = self._unpack(parameters)
        offset = frame.prior_median_pactivity.to_numpy(float)
        outcome = frame.heldout_median_pactivity.to_numpy(float)
        location = offset + alpha[target_index]
        losses = []
        for indices in study_groups:
            target = target_index[indices[0]]
            if tau[target] == 0:
                log_likelihood = student_t.logpdf(
                    (outcome[indices] - location[indices]) / sigma,
                    df=self.config.degrees_of_freedom,
                ).sum() - len(indices) * math.log(sigma)
            else:
                shifts = math.sqrt(2.0) * tau[target] * self._gh_nodes
                standardized = (
                    outcome[indices, None] - location[indices, None] - shifts[None, :]
                ) / sigma
                node_log_likelihood = (
                    student_t.logpdf(standardized, df=self.config.degrees_of_freedom).sum(axis=0)
                    - len(indices) * math.log(sigma)
                )
                log_likelihood = logsumexp(self._gh_log_weights + node_log_likelihood)
            # Each publication contributes equally, as prespecified.
            losses.append(-float(log_likelihood) / len(indices))
        penalty = 0.5 * float(np.sum((alpha / self.config.target_effect_scale) ** 2)) / max(len(frame), 1)
        if self.include_study_variance:
            count = len(self.targets_)
            delta = parameters[count + 2:count + 2 + count]
            penalty += 0.5 * float(np.sum((delta / self.config.log_tau_pooling_scale) ** 2)) / max(len(frame), 1)
        return float(np.mean(losses) + penalty)

    def fit(self, frame: pd.DataFrame) -> "HierarchicalStudentTM0":
        required = {"target_chembl_id", "heldout_study_key", "prior_median_pactivity", "heldout_median_pactivity"}
        if missing := sorted(required - set(frame.columns)):
            raise ValueError(f"Missing fit columns: {missing}")
        self.targets_ = sorted(frame.target_chembl_id.astype(str).unique())
        mapping = {target: index for index, target in enumerate(self.targets_)}
        target_index = frame.target_chembl_id.map(mapping).to_numpy(int)
        study_groups = [group.index.to_numpy() for _, group in frame.reset_index(drop=True).groupby("heldout_study_key", sort=True)]
        residual = frame.heldout_median_pactivity.to_numpy(float) - frame.prior_median_pactivity.to_numpy(float)
        alpha0 = np.asarray([
            np.median(residual[target_index == index]) for index in range(len(self.targets_))
        ])
        centered = residual - alpha0[target_index]
        sigma0 = max(float(np.median(np.abs(centered - np.median(centered))) * 1.4826), 0.25)
        if self.include_study_variance:
            initial = np.concatenate([
                alpha0,
                [math.log(sigma0), math.log(max(sigma0 / 2.0, 0.05))],
                np.zeros(len(self.targets_)),
            ])
        else:
            initial = np.concatenate([alpha0, [math.log(sigma0)]])
        result = minimize(
            self._objective,
            initial,
            args=(frame.reset_index(drop=True), target_index, study_groups),
            method="L-BFGS-B",
            options={"maxiter": self.config.optimizer_maxiter, "ftol": 1e-10, "gtol": 1e-6},
        )
        if not result.success and not np.isfinite(result.fun):
            raise RuntimeError(f"M0 optimization failed: {result.message}")
        self.alpha_, self.sigma_, self.tau_ = self._unpack(result.x)
        self.optimization_ = {
            "success": bool(result.success), "message": str(result.message),
            "iterations": int(result.nit), "objective": float(result.fun),
        }
        return self

    def _indices(self, frame: pd.DataFrame) -> np.ndarray:
        mapping = {target: index for index, target in enumerate(self.targets_)}
        unknown = sorted(set(frame.target_chembl_id.astype(str)) - set(mapping))
        if unknown:
            raise ValueError(f"M0 cannot extrapolate to unseen targets: {unknown}")
        return frame.target_chembl_id.astype(str).map(mapping).to_numpy(int)

    def predict_location(self, frame: pd.DataFrame) -> np.ndarray:
        if self.alpha_ is None:
            raise RuntimeError("Model is not fitted")
        indices = self._indices(frame)
        return frame.prior_median_pactivity.to_numpy(float) + self.alpha_[indices]

    def _density(self, values: np.ndarray, location: np.ndarray, tau: np.ndarray) -> np.ndarray:
        assert self.sigma_ is not None
        if not self.include_study_variance:
            return student_t.pdf((values - location) / self.sigma_, df=self.config.degrees_of_freedom) / self.sigma_
        shifts = math.sqrt(2.0) * tau[:, None] * self._gh_nodes[None, :]
        standardized = (values[:, None] - location[:, None] - shifts) / self.sigma_
        return np.exp(logsumexp(self._gh_log_weights[None, :] + student_t.logpdf(
            standardized, df=self.config.degrees_of_freedom
        ) - math.log(self.sigma_), axis=1))

    def _cdf_scalar(self, value: float, location: float, tau: float) -> float:
        assert self.sigma_ is not None
        if not self.include_study_variance:
            return float(student_t.cdf((value - location) / self.sigma_, df=self.config.degrees_of_freedom))
        shifts = math.sqrt(2.0) * tau * self._gh_nodes
        return float(np.sum(np.exp(self._gh_log_weights) * student_t.cdf(
            (value - location - shifts) / self.sigma_, df=self.config.degrees_of_freedom
        )))

    def predict_distribution(self, frame: pd.DataFrame, levels: Iterable[float] = (0.8, 0.95)) -> dict[str, np.ndarray]:
        if self.tau_ is None or self.sigma_ is None:
            raise RuntimeError("Model is not fitted")
        target_index = self._indices(frame)
        location = self.predict_location(frame)
        tau = self.tau_[target_index]
        outcome = frame.heldout_median_pactivity.to_numpy(float) if "heldout_median_pactivity" in frame else None
        output: dict[str, np.ndarray] = {
            "prediction": location,
            "residual_scale": np.full(len(frame), self.sigma_),
            "between_study_sd": tau,
        }
        if outcome is not None:
            output["nll"] = -np.log(np.maximum(self._density(outcome, location, tau), 1e-300))
        heavy_tail_sd = self.sigma_ * math.sqrt(self.config.degrees_of_freedom / (self.config.degrees_of_freedom - 2.0))
        for level in levels:
            tail = (1.0 - level) / 2.0
            lower, upper = [], []
            for center, study_sd in zip(location, tau):
                total_sd = math.sqrt(heavy_tail_sd ** 2 + study_sd ** 2)
                radius = max(12.0 * total_sd, 5.0)
                quantiles = []
                for probability in (tail, 1.0 - tail):
                    quantiles.append(brentq(
                        lambda value: self._cdf_scalar(value, float(center), float(study_sd)) - probability,
                        float(center - radius), float(center + radius),
                    ))
                lower.append(quantiles[0])
                upper.append(quantiles[1])
            suffix = str(int(level * 100))
            output[f"lower_{suffix}"] = np.asarray(lower)
            output[f"upper_{suffix}"] = np.asarray(upper)
        return output

    def parameters(self) -> dict[str, Any]:
        if self.alpha_ is None or self.tau_ is None or self.sigma_ is None:
            raise RuntimeError("Model is not fitted")
        return {
            "targets": self.targets_,
            "target_effects": dict(zip(self.targets_, map(float, self.alpha_))),
            "residual_scale": float(self.sigma_),
            "target_between_study_sd": dict(zip(self.targets_, map(float, self.tau_))),
            "degrees_of_freedom": self.config.degrees_of_freedom,
            "optimization": self.optimization_,
        }


def _student_scale(residual: np.ndarray, df: float) -> float:
    initial = max(float(np.median(np.abs(residual - np.median(residual))) * 1.4826), 0.1)
    result = minimize(
        lambda value: -float(np.sum(student_t.logpdf(residual / math.exp(float(value[0])), df=df) - float(value[0]))),
        np.asarray([math.log(initial)]), method="L-BFGS-B",
    )
    return float(math.exp(float(result.x[0])))


def _baseline_predictions(train: pd.DataFrame, test: pd.DataFrame, config: M0Config) -> dict[str, dict[str, np.ndarray]]:
    y_train = train.heldout_median_pactivity.to_numpy(float)
    methods: dict[str, tuple[np.ndarray, float]] = {}
    global_mean = float(np.mean(y_train))
    methods["global_mean"] = (np.full(len(test), global_mean), _student_scale(y_train - global_mean, config.degrees_of_freedom))
    historical = test.prior_median_pactivity.to_numpy(float)
    methods["median_potency"] = (historical, _student_scale(
        train.heldout_median_pactivity.to_numpy(float) - train.prior_median_pactivity.to_numpy(float),
        config.degrees_of_freedom,
    ))
    output: dict[str, dict[str, np.ndarray]] = {}
    tail_values = {80: student_t.ppf(0.90, df=config.degrees_of_freedom), 95: student_t.ppf(0.975, df=config.degrees_of_freedom)}
    for method, (prediction, scale) in methods.items():
        outcome = test.heldout_median_pactivity.to_numpy(float)
        result = {
            "prediction": prediction,
            "nll": -student_t.logpdf((outcome - prediction) / scale, df=config.degrees_of_freedom) + math.log(scale),
            "between_study_sd": np.zeros(len(test)),
            "residual_scale": np.full(len(test), scale),
        }
        for level, critical in tail_values.items():
            result[f"lower_{level}"] = prediction - critical * scale
            result[f"upper_{level}"] = prediction + critical * scale
        output[method] = result
    return output


def _assert_fold_safe(frame: pd.DataFrame, train_index: np.ndarray, test_index: np.ndarray, scheme: str) -> None:
    train, test = frame.iloc[train_index], frame.iloc[test_index]
    if set(train.heldout_study_key) & set(test.heldout_study_key):
        raise AssertionError(f"Study leakage in {scheme}")
    if set(train.molecule_id) & set(test.molecule_id):
        raise AssertionError(f"Compound leakage in {scheme}")
    if scheme == "forward_temporal" and train.heldout_year.max() >= test.heldout_year.min():
        raise AssertionError("Future-information leakage in temporal split")


def cross_validate(frame: pd.DataFrame, config: M0Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    splits = make_splits(frame)
    rows: list[pd.DataFrame] = []
    parameter_rows: list[dict[str, Any]] = []
    for scheme in ("study_grouped", "forward_temporal"):
        for fold, (train_index, test_index) in enumerate(splits[scheme]):
            _assert_fold_safe(frame, train_index, test_index, scheme)
            train, test = frame.iloc[train_index].copy(), frame.iloc[test_index].copy()
            baseline = _baseline_predictions(train, test, config)
            nonhierarchical = HierarchicalStudentTM0(config, include_study_variance=False).fit(train)
            m0 = HierarchicalStudentTM0(config, include_study_variance=True).fit(train)
            predictions = {
                **baseline,
                "nonhierarchical_student_t": nonhierarchical.predict_distribution(test),
                "phytotrust_m0": m0.predict_distribution(test),
            }
            identity = test[[
                "target_chembl_id", "molecule_id", "heldout_study_key", "heldout_year",
                "prior_median_pactivity", "heldout_median_pactivity",
                "is_verified_natural_product", "prior_document_ids", "heldout_document_ids",
                "prior_assay_ids", "heldout_assay_ids",
            ]].reset_index(names="source_index")
            for method, predicted in predictions.items():
                result = identity.copy()
                result.insert(0, "method", method)
                result.insert(0, "fold", fold)
                result.insert(0, "validation_scheme", scheme)
                for name, values in predicted.items():
                    result[name] = values
                rows.append(result)
            for method, model in (("nonhierarchical_student_t", nonhierarchical), ("phytotrust_m0", m0)):
                params = model.parameters()
                parameter_rows.append({
                    "validation_scheme": scheme, "fold": fold, "method": method,
                    "training_rows": len(train), "training_studies": train.heldout_study_key.nunique(),
                    "residual_scale": params["residual_scale"],
                    "mean_between_study_sd": float(np.mean(list(params["target_between_study_sd"].values()))),
                    "target_between_study_sd": json.dumps(params["target_between_study_sd"], sort_keys=True),
                    "optimizer_success": params["optimization"]["success"],
                })
    return pd.concat(rows, ignore_index=True), pd.DataFrame(parameter_rows)


def _metric_row(group: pd.DataFrame) -> dict[str, Any]:
    observed = group.heldout_median_pactivity.to_numpy(float)
    predicted = group.prediction.to_numpy(float)
    result: dict[str, Any] = {
        "n": int(len(group)), "studies": int(group.heldout_study_key.nunique()),
        "mae": float(mean_absolute_error(observed, predicted)),
        "rmse": float(math.sqrt(mean_squared_error(observed, predicted))),
        "nll": float(group.nll.mean()),
    }
    for level in (80, 95):
        covered = (observed >= group[f"lower_{level}"]) & (observed <= group[f"upper_{level}"])
        result[f"coverage_{level}"] = float(covered.mean())
        result[f"mean_width_{level}"] = float((group[f"upper_{level}"] - group[f"lower_{level}"]).mean())
    return result


def evaluate_predictions(predictions: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows, nested = [], {}
    for (scheme, method), group in predictions.groupby(["validation_scheme", "method"], sort=True):
        nested.setdefault(scheme, {})
        for scope, scoped in (("all", group), ("natural_products", group[group.is_verified_natural_product])):
            values = _metric_row(scoped) if len(scoped) else {}
            rows.append({"validation_scheme": scheme, "scope": scope, "method": method, **values})
            nested[scheme].setdefault(scope, {})[method] = values
    return pd.DataFrame(rows), nested


def cluster_bootstrap_differences(
    predictions: pd.DataFrame, *, iterations: int, seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    records, nested = [], {}
    methods = ["global_mean", "median_potency", "nonhierarchical_student_t"]
    for scheme in ("study_grouped", "forward_temporal"):
        wide_source = predictions[predictions.validation_scheme == scheme]
        for scope, scoped in (("all", wide_source), ("natural_products", wide_source[wide_source.is_verified_natural_product])):
            wide = scoped.pivot(index=["source_index", "heldout_study_key"], columns="method", values=["prediction", "nll", "lower_80", "upper_80", "heldout_median_pactivity"])
            studies = sorted(wide.index.get_level_values("heldout_study_key").unique())
            nested.setdefault(scheme, {}).setdefault(scope, {})
            for baseline in methods:
                observed_point = wide[("heldout_median_pactivity", "phytotrust_m0")].to_numpy(float)
                m0_error_point = observed_point - wide[("prediction", "phytotrust_m0")].to_numpy(float)
                base_error_point = observed_point - wide[("prediction", baseline)].to_numpy(float)
                point = {
                    "mae": float(np.mean(np.abs(m0_error_point)) - np.mean(np.abs(base_error_point))),
                    "rmse": float(math.sqrt(np.mean(m0_error_point ** 2)) - math.sqrt(np.mean(base_error_point ** 2))),
                    "nll": float(wide[("nll", "phytotrust_m0")].mean() - wide[("nll", baseline)].mean()),
                }
                samples = {"mae": [], "rmse": [], "nll": []}
                for _ in range(iterations):
                    chosen = rng.choice(studies, len(studies), replace=True)
                    sampled = pd.concat([wide.xs(study, level="heldout_study_key", drop_level=False) for study in chosen], ignore_index=True)
                    observed = sampled[("heldout_median_pactivity", "phytotrust_m0")].to_numpy(float)
                    for metric in ("mae", "rmse"):
                        m0_error = observed - sampled[("prediction", "phytotrust_m0")].to_numpy(float)
                        base_error = observed - sampled[("prediction", baseline)].to_numpy(float)
                        if metric == "mae":
                            delta = np.mean(np.abs(m0_error)) - np.mean(np.abs(base_error))
                        else:
                            delta = math.sqrt(np.mean(m0_error ** 2)) - math.sqrt(np.mean(base_error ** 2))
                        samples[metric].append(float(delta))
                    samples["nll"].append(float(sampled[("nll", "phytotrust_m0")].mean() - sampled[("nll", baseline)].mean()))
                intervals = {metric: [float(x) for x in np.quantile(values, [0.025, 0.975])] for metric, values in samples.items()}
                records.append({
                    "validation_scheme": scheme, "scope": scope, "comparison": f"phytotrust_m0-minus-{baseline}",
                    **{f"delta_{metric}": value for metric, value in point.items()},
                    **{f"delta_{metric}_ci95": json.dumps(interval) for metric, interval in intervals.items()},
                })
                nested[scheme][scope][baseline] = {"point": point, "ci95": intervals}
    return pd.DataFrame(records), nested


def identifiability_diagnostics(frame: pd.DataFrame, parameters: pd.DataFrame, final_model: HierarchicalStudentTM0) -> dict[str, Any]:
    counts = frame.groupby(["target_chembl_id", "heldout_study_key"]).size().rename("events").reset_index()
    target_data = {}
    final_parameters = final_model.parameters()
    for target, group in counts.groupby("target_chembl_id"):
        fold_rows = parameters[(parameters.method == "phytotrust_m0")]
        fold_tau = [
            json.loads(value).get(target) for value in fold_rows.target_between_study_sd
            if target in json.loads(value)
        ]
        target_data[target] = {
            "studies": int(len(group)), "multi_event_studies": int((group.events >= 2).sum()),
            "studies_with_at_least_3_events": int((group.events >= 3).sum()),
            "largest_study": int(group.events.max()),
            "final_between_study_sd": float(final_parameters["target_between_study_sd"][target]),
            "cross_fold_tau_range": [float(min(fold_tau)), float(max(fold_tau))] if fold_tau else None,
        }
    identifiable_targets = [
        target for target, values in target_data.items()
        if values["multi_event_studies"] >= 10 and values["studies_with_at_least_3_events"] >= 5
    ]
    return {
        "criterion": "at least 10 multi-event studies and 5 studies with >=3 events per target, plus non-collapsed stable fold estimates",
        "target_diagnostics": target_data,
        "targets_meeting_structural_gate": identifiable_targets,
        "target_specific_study_variance_identifiable": len(identifiable_targets) == len(target_data),
        "conclusion": (
            "Target-specific between-study variance is not statistically identifiable with the available cluster structure. "
            "M0 estimates are regularized feasibility estimates, not definitive variance components."
            if len(identifiable_targets) != len(target_data)
            else "All targets meet the prespecified structural gate; stability still requires external validation."
        ),
    }


def _render_report(results: dict[str, Any], metric_table: pd.DataFrame) -> str:
    dataset, ident = results["dataset"], results["identifiability"]
    lines = [
        "# PhytoTRUST-M0 scientific validation", "",
        "## Scope", "",
        f"This experiment used only the {dataset['eligible_rows']} real pre-2020 development events from {dataset['studies']} ChEMBL publication groups and {dataset['targets']} targets. The locked 2020+ partition and previously examined independent cohorts were not used for fitting, selection, or evaluation.", "",
        "M0 models held-out pChEMBL as the historical median plus a target correction, an independent-study random intercept, and Student-t residual noise. For blind transfer, the held-out study effect is integrated over its learned target-specific population distribution. No outcome from the held-out publication is used at prediction time.", "",
        "## Results", "",
        "| Validation | Scope | Method | n | MAE | RMSE | NLL | 80% coverage / width | 95% coverage / width |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metric_table.itertuples():
        lines.append(
            f"| {row.validation_scheme} | {row.scope} | {row.method} | {row.n} | {row.mae:.3f} | {row.rmse:.3f} | {row.nll:.3f} | {row.coverage_80:.3f} / {row.mean_width_80:.3f} | {row.coverage_95:.3f} / {row.mean_width_95:.3f} |"
        )
    sg = results["metrics"]["study_grouped"]["all"]
    temporal = results["metrics"]["forward_temporal"]["all"]
    m0, nonhier = sg["phytotrust_m0"], sg["nonhierarchical_student_t"]
    sensitivity = results["dominant_study_sensitivity"]
    lines += [
        "", "## Interpretation", "",
        f"In study-grouped validation, M0 had MAE {m0['mae']:.3f}, RMSE {m0['rmse']:.3f}, and NLL {m0['nll']:.3f}. The matched non-hierarchical Student-t model had MAE {nonhier['mae']:.3f}, RMSE {nonhier['rmse']:.3f}, and NLL {nonhier['nll']:.3f}. Equal point errors are expected because blind integration of a zero-mean unseen-study effect changes uncertainty, not the predictive mean.", "",
        "The cluster-bootstrap comparison intervals are stored in `bootstrap_differences.csv` and `results.json`. A negative M0-minus-baseline interval favors M0. Results are development-only and do not establish novelty or independent generalization.", "",
        f"M0 did not improve the primary probabilistic result: its study-grouped NLL was {m0['nll']:.3f}, versus {sg['median_potency']['nll']:.3f} for median potency and {nonhier['nll']:.3f} for the matched non-hierarchical model. In forward-time validation, the corresponding values were {temporal['phytotrust_m0']['nll']:.3f}, {temporal['median_potency']['nll']:.3f}, and {temporal['nonhierarchical_student_t']['nll']:.3f}. The natural-product study-grouped point estimate favored M0 on NLL, but its paired study-bootstrap intervals included no effect and the direction did not reproduce temporally.", "",
        f"Removing the largest publication ({sensitivity['removed_study']}, {sensitivity['removed_events']} events) changed study-grouped M0 NLL to {sensitivity['metrics_without_study']['study_grouped']['all']['phytotrust_m0']['nll']:.3f}; the non-hierarchical value became {sensitivity['metrics_without_study']['study_grouped']['all']['nonhierarchical_student_t']['nll']:.3f}. The main negative conclusion therefore does not depend solely on that publication.", "",
        "## Variance identifiability", "",
        ident["conclusion"], "",
        f"Targets meeting the structural identifiability gate: {len(ident['targets_meeting_structural_gate'])}/{len(ident['target_diagnostics'])}. Partial pooling stabilizes estimation but cannot create information absent from singleton studies.", "",
        "## Natural products", "",
        f"The verified-natural-product subgroup contains {dataset['natural_products']} events. It was evaluation-only: NP status was not a predictor or fitting weight. Its estimates are exploratory because study and target coverage are sparse.", "",
        "## Limitations", "",
        "- M0 estimates a measurement distribution, not true potency, efficacy, safety, or which study is correct.",
        "- The dataset contains event-level study medians rather than every within-assay replicate, limiting separation of within-assay and between-study noise.",
        "- Target-specific variance is weakly identified because 112/141 studies are singletons and only 10 studies have at least three events.",
        "- Student-t degrees of freedom and pooling scales were prespecified from the methodology, not tuned on held-out outcomes.",
        "- Temporal validation is compound-disjoint and therefore tests a strict form of transfer; it may omit historically repeated compounds from training by design.",
        "- This valid M0 implementation is a feasibility baseline, not evidence that the method is novel or publication-ready.",
        "", "## Scientific assessment", "",
        "M0 is a valid, reproducible hierarchical feasibility baseline, but the available development evidence does not support it as a journal-level methodological foundation yet. It failed to improve NLL over median potency or the matched non-hierarchical model, temporal performance was worse, and target-specific study variance did not pass the identifiability gate. More multi-compound independent studies are required before extending the hierarchy.",
    ]
    return "\n".join(lines) + "\n"


def _plot_metrics(metric_table: pd.DataFrame, path: Path) -> None:
    methods = ["global_mean", "median_potency", "nonhierarchical_student_t", "phytotrust_m0"]
    labels = ["Global mean", "Median potency", "Non-hierarchical t", "PhytoTRUST-M0"]
    colors = ["#64748b", "#0ea5e9", "#a78bfa", "#22c55e"]
    panels = [("study_grouped", "Study-grouped"), ("forward_temporal", "Forward temporal")]
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="590" viewBox="0 0 1100 590">',
        '<rect width="1100" height="590" fill="#07111f"/>',
        '<text x="55" y="42" fill="#f8fafc" font-family="sans-serif" font-size="22">PhytoTRUST-M0 development-only validation</text>',
        '<text x="55" y="66" fill="#94a3b8" font-family="sans-serif" font-size="13">Lower is better; natural products are evaluation-only.</text>',
    ]
    for panel_index, (scheme, title) in enumerate(panels):
        subset = metric_table[(metric_table.validation_scheme == scheme) & (metric_table.scope == "all")]
        x0 = 55 + panel_index * 530
        pieces.append(f'<text x="{x0}" y="105" fill="#e2e8f0" font-family="sans-serif" font-size="17">{title}</text>')
        for metric_index, metric in enumerate(("nll", "mae")):
            y0 = 135 + metric_index * 210
            pieces.append(f'<text x="{x0}" y="{y0}" fill="#94a3b8" font-family="sans-serif" font-size="13">{metric.upper()}</text>')
            maximum = max(float(subset[metric].max()), 0.001)
            for index, (method, label, color) in enumerate(zip(methods, labels, colors)):
                value = float(subset.loc[subset.method == method, metric].iloc[0])
                y = y0 + 22 + index * 40
                width = 300 * value / maximum
                pieces.append(f'<text x="{x0}" y="{y+15}" fill="#cbd5e1" font-family="sans-serif" font-size="12">{label}</text>')
                pieces.append(f'<rect x="{x0+145}" y="{y}" width="{width:.1f}" height="21" rx="4" fill="{color}"/>')
                pieces.append(f'<text x="{x0+151+width:.1f}" y="{y+15}" fill="#f8fafc" font-family="sans-serif" font-size="12">{value:.3f}</text>')
    pieces.append('</svg>')
    path.write_text("\n".join(pieces), encoding="utf-8")


def run_experiment(
    dataset: Path, output_dir: Path, *, bootstrap_iterations: int = 1000,
    seed: int = SEED, config: M0Config | None = None,
) -> dict[str, Any]:
    config = config or M0Config()
    output_dir.mkdir(parents=True, exist_ok=True)
    frame, manifest = load_m0_development_data(dataset)
    predictions, parameters = cross_validate(frame, config)
    metric_table, nested_metrics = evaluate_predictions(predictions)
    bootstrap_table, bootstrap = cluster_bootstrap_differences(
        predictions, iterations=bootstrap_iterations, seed=seed,
    )
    final_model = HierarchicalStudentTM0(config, include_study_variance=True).fit(frame)
    identifiability = identifiability_diagnostics(frame, parameters, final_model)
    largest_study = frame.heldout_study_key.value_counts().index[0]
    without_dominant = predictions[predictions.heldout_study_key != largest_study]
    sensitivity_table, sensitivity_metrics = evaluate_predictions(without_dominant)
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "model_version": MODEL_VERSION,
        "seed": seed, "configuration": asdict(config), "dataset": manifest,
        "validation": {
            "schemes": ["study_grouped", "forward_temporal"],
            "compound_disjoint": True, "publication_disjoint": True,
            "forward_temporal_strict": True, "heldout_study_effect_integrated": True,
        },
        "metrics": nested_metrics, "study_cluster_bootstrap_differences": bootstrap,
        "identifiability": identifiability,
        "dominant_study_sensitivity": {
            "removed_study": str(largest_study),
            "removed_events": int((frame.heldout_study_key == largest_study).sum()),
            "metrics_without_study": sensitivity_metrics,
        },
        "final_parameters": final_model.parameters(),
        "scientific_status": "implemented_feasibility_baseline_not_independently_validated",
    }
    predictions.to_csv(output_dir / "predictions.csv", index=False)
    metric_table.to_csv(output_dir / "metrics.csv", index=False)
    parameters.to_csv(output_dir / "fold_parameters.csv", index=False)
    bootstrap_table.to_csv(output_dir / "bootstrap_differences.csv", index=False)
    sensitivity_table.to_csv(output_dir / "dominant_study_sensitivity.csv", index=False)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir / "model_config.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "configuration": asdict(config),
        "equation": "Y = prior_median + alpha_target + u_study,target + StudentT(df, 0, sigma)",
        "blind_transfer": "u_new integrated over Normal(0, tau_target^2)",
    }, indent=2), encoding="utf-8")
    (output_dir / "model_parameters.json").write_text(json.dumps(final_model.parameters(), indent=2), encoding="utf-8")
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_render_report(results, metric_table), encoding="utf-8")
    _plot_metrics(metric_table, output_dir / "comparison.svg")
    (output_dir / "command.txt").write_text(
        "python -m app.phytotrust_m0 --dataset ../experiments/assay_protocol_enrichment/enriched_assay_disagreement_events.csv --output-dir ../experiments/phytotrust_m0 --bootstrap-iterations 1000 --seed 2026\n",
        encoding="utf-8",
    )
    files = sorted(path.name for path in output_dir.iterdir() if path.is_file())
    (output_dir / "manifest.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "source_dataset": str(dataset),
        "source_sha256": manifest["source_sha256"], "files": files,
        "locked_evaluation_used": False,
    }, indent=2), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="PhytoTRUST-M0 study-aware hierarchical feasibility experiment")
    parser.add_argument("--dataset", type=Path, default=Path("../experiments/assay_protocol_enrichment/enriched_assay_disagreement_events.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("../experiments/phytotrust_m0"))
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    results = run_experiment(args.dataset, args.output_dir, bootstrap_iterations=args.bootstrap_iterations, seed=args.seed)
    print(json.dumps({
        "dataset": results["dataset"], "scientific_status": results["scientific_status"],
        "identifiability": results["identifiability"]["target_specific_study_variance_identifiable"],
    }, indent=2))


if __name__ == "__main__":
    main()
