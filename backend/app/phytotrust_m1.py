from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp
from scipy.stats import t as student_t
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

from .assay_protocol import EXTRACTION_VERSION, SOURCE_BASE, SOURCE_NAME, extract_protocol_record, load_and_extract
from .chembl import iter_assays_by_ids
from .phytotrust_m0 import HierarchicalStudentTM0, M0Config, _student_scale
from .phytotrust_m0_connected import (
    EXPECTED_INPUT_NAME,
    aggregate_compound_study,
    assert_fold_provenance,
    load_connected_development,
)


SEED = 2026
MODEL_VERSION = "phytotrust-m1-v1"
MEAN_FEATURES = (
    "history_mean_minus_median", "history_mad", "log1p_prior_studies",
    "molecular_weight_scaled", "logp_scaled",
)
SCALE_FEATURES = (
    "history_mad", "inverse_sqrt_prior_studies", "protocol_missing_fraction",
    "protocol_category_diversity", "history_year_span_scaled",
)
PROTOCOL_SCALE_FEATURES = {"protocol_missing_fraction", "protocol_category_diversity"}
ALL_TABULAR_FEATURES = tuple(dict.fromkeys((*MEAN_FEATURES, *SCALE_FEATURES)))
METHODS = (
    "median_potency", "student_t_regression", "xgboost", "phytotrust_m0",
    "m1_no_protocol", "m1_homoscedastic", "m1_no_study_variance", "phytotrust_m1",
)


@dataclass(frozen=True)
class M1Config:
    degrees_of_freedom: float = 4.0
    quadrature_nodes: int = 15
    mean_coefficient_scale: float = 0.35
    log_scale_coefficient_scale: float = 0.5
    minimum_scale: float = 1e-3
    optimizer_maxiter: int = 900
    optimizer_maxfun: int = 100000


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _joined(values: Iterable[Any]) -> str:
    return ";".join(sorted({str(value) for value in values if pd.notna(value) and str(value)}))


async def fetch_development_assays(
    dataset: Path, raw_output: Path, feature_output: Path, *, batch_size: int = 50, max_retries: int = 5,
) -> dict[str, Any]:
    frame, development_manifest = load_connected_development(dataset)
    assay_ids = sorted(frame.assay_id.astype(str).unique())
    raw_output.parent.mkdir(parents=True, exist_ok=True)
    part = raw_output.with_suffix(raw_output.suffix + ".part")
    retrieved_at = datetime.now(UTC).isoformat()
    returned: set[str] = set()
    digest = hashlib.sha256()
    with part.open("w", encoding="utf-8", newline="\n") as handle:
        async for record, request_url in iter_assays_by_ids(
            assay_ids, batch_size=batch_size, max_retries=max_retries,
        ):
            identifier = str(record.get("assay_chembl_id") or "")
            if not identifier or identifier in returned:
                continue
            raw_hash = hashlib.sha256(json.dumps(
                record, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
            ).encode("utf-8")).hexdigest()
            envelope = {
                "source": SOURCE_NAME, "source_endpoint": request_url,
                "retrieved_at": retrieved_at, "assay_chembl_id": identifier,
                "raw_record_sha256": raw_hash, "record": record,
            }
            line = json.dumps(envelope, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
            handle.write(line)
            digest.update(line.encode("utf-8"))
            returned.add(identifier)
    part.replace(raw_output)
    missing = sorted(set(assay_ids) - returned)
    if missing:
        raise RuntimeError(f"ChEMBL did not return {len(missing)} development assays: {missing[:5]}")
    features = load_and_extract(raw_output)
    features.to_csv(feature_output, index=False)
    manifest = {
        "source": SOURCE_NAME, "source_base": SOURCE_BASE, "retrieved_at": retrieved_at,
        "development_dataset": str(dataset.resolve()),
        "development_dataset_sha256": development_manifest["source_sha256"],
        "requested_assays": len(assay_ids), "returned_assays": len(returned),
        "missing_assay_ids": missing, "batch_size": batch_size, "max_retries": max_retries,
        "raw_output": raw_output.name, "raw_sha256": digest.hexdigest(),
        "feature_output": feature_output.name, "feature_sha256": _sha256(feature_output),
        "extraction_version": EXTRACTION_VERSION,
        "quarantined_files_opened": False,
    }
    raw_output.with_suffix(raw_output.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8",
    )
    return manifest


def load_protocol_features(path: Path, expected_assays: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    lowered = str(path.resolve()).lower()
    if "quarantined" in lowered or "future_evaluation" in lowered:
        raise ValueError("Evaluation protocol metadata is forbidden")
    frame = pd.read_csv(path)
    if frame.assay_chembl_id.astype(str).duplicated().any():
        raise AssertionError("Duplicate protocol metadata rows")
    found = set(frame.assay_chembl_id.astype(str))
    missing = sorted(expected_assays - found)
    extra = sorted(found - expected_assays)
    if missing or extra:
        raise AssertionError(f"Protocol assay mismatch: missing={missing[:5]}, extra={extra[:5]}")
    return frame, {
        "path": str(path.resolve()), "sha256": _sha256(path), "assays": len(frame),
        "missing_assays": missing, "extra_assays": extra,
    }


def _tokens(value: Any) -> set[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    return {item.strip() for item in str(value).split(";") if item.strip()}


def _molecular_descriptors(smiles: str) -> dict[str, float]:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid standardized SMILES in curated development dataset")
    return {
        "molecular_weight_scaled": float(Descriptors.MolWt(molecule) / 500.0),
        "logp_scaled": float(Crippen.MolLogP(molecule) / 5.0),
    }


def enrich_edges(measurements: pd.DataFrame, protocol: pd.DataFrame) -> pd.DataFrame:
    base = aggregate_compound_study(measurements)
    molecule_smiles = measurements.groupby("molecule_id").standardized_smiles.first()
    protocol_index = protocol.set_index(protocol.assay_chembl_id.astype(str), drop=False)
    missing_fields = (
        "assay_cell_type", "assay_tissue", "assay_subcellular_fraction", "assay_strain",
        "assay_parameter_types", "variant_mutation", "detection_categories",
        "methodology_categories", "temperature_mentions", "ph_mentions", "duration_mentions",
        "substrate_mentions_json",
    )
    token_fields = (
        "methodology_categories", "detection_categories", "biological_system_text_categories",
        "assay_cell_type", "assay_tissue", "assay_subcellular_fraction", "assay_parameter_types",
        "temperature_mentions", "ph_mentions", "duration_mentions",
    )
    records = []
    descriptor_cache: dict[str, dict[str, float]] = {}
    for row in base.itertuples(index=False):
        assay_ids = [value for value in str(row.assay_ids).split(";") if value]
        assays = protocol_index.loc[assay_ids]
        if isinstance(assays, pd.Series):
            assays = assays.to_frame().T
        missing_cells = 0
        total_cells = len(assays) * len(missing_fields)
        for field in missing_fields:
            values = assays[field] if field in assays else pd.Series([None] * len(assays))
            missing_cells += int(sum(value is None or pd.isna(value) or str(value).strip() in ("", "[]") for value in values))
        tokens: set[str] = set()
        for field in token_fields:
            if field not in assays:
                continue
            for value in assays[field]:
                for token in _tokens(value):
                    tokens.add(f"{field}:{token}")
        molecule_id = str(row.molecule_id)
        if molecule_id not in descriptor_cache:
            descriptor_cache[molecule_id] = _molecular_descriptors(str(molecule_smiles.loc[molecule_id]))
        record = row._asdict()
        record.update(descriptor_cache[molecule_id])
        record.update({
            "protocol_assay_count": len(assay_ids),
            "protocol_missing_fraction": float(missing_cells / total_cells) if total_cells else 1.0,
            "protocol_category_count": len(tokens),
            "protocol_tokens": _joined(tokens),
            "protocol_direct_completeness": float(pd.to_numeric(
                assays.direct_field_completeness, errors="coerce",
            ).mean()),
            "protocol_source_assay_ids": _joined(assay_ids),
        })
        records.append(record)
    return pd.DataFrame(records).sort_values(
        ["document_year", "document_id", "molecule_id"],
    ).reset_index(drop=True)


def _history_features(prior: pd.DataFrame) -> dict[str, float]:
    values = prior.study_median_pactivity.to_numpy(float)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median))) if len(values) > 1 else 0.0
    tokens: set[str] = set()
    for value in prior.protocol_tokens:
        tokens.update(_tokens(value))
    assay_count = max(int(prior.protocol_assay_count.sum()), 1)
    return {
        "prior_median_pactivity": median,
        "history_mean_minus_median": float(np.mean(values) - median),
        "history_mad": mad,
        "log1p_prior_studies": float(math.log1p(prior.document_id.nunique())),
        "inverse_sqrt_prior_studies": float(1.0 / math.sqrt(prior.document_id.nunique())),
        "protocol_missing_fraction": float(np.average(
            prior.protocol_missing_fraction, weights=prior.protocol_assay_count,
        )),
        "protocol_category_diversity": float(len(tokens) / assay_count),
        "history_year_span_scaled": float((prior.document_year.max() - prior.document_year.min()) / 10.0),
    }


def construct_feature_events(
    source_edges: pd.DataFrame, heldout_edges: pd.DataFrame | None = None, *, history_mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if heldout_edges is None:
        candidates = source_edges
        grouped = {key: value for key, value in source_edges.groupby("molecule_id", sort=False)}
        leave_one_out = True
    else:
        candidates = heldout_edges
        grouped = {key: value for key, value in source_edges.groupby("molecule_id", sort=False)}
        leave_one_out = False
    eligible_candidates = candidates[candidates.molecule_id.isin(grouped)]
    bridge_counts = eligible_candidates.groupby("document_id").molecule_id.nunique()
    for heldout in candidates.itertuples(index=False):
        compound = grouped.get(heldout.molecule_id)
        if compound is None:
            continue
        prior = compound[compound.document_id != heldout.document_id] if leave_one_out else compound
        if prior.empty:
            continue
        feature_values = _history_features(prior)
        prior_documents = _joined(prior.document_id)
        prior_assays = _joined(";".join(prior.assay_ids).split(";"))
        prior_activities = _joined(";".join(prior.activity_ids).split(";"))
        if str(heldout.document_id) in prior_documents.split(";"):
            raise AssertionError("Held-out publication entered M1 history")
        if set(str(heldout.activity_ids).split(";")) & set(prior_activities.split(";")):
            raise AssertionError("Held-out activity entered M1 history")
        count = int(prior.document_id.nunique())
        row = {
            "target_chembl_id": heldout.target_chembl_id,
            "molecule_id": heldout.molecule_id,
            "heldout_study_key": heldout.document_id,
            "heldout_year": int(heldout.document_year),
            "heldout_median_pactivity": float(heldout.study_median_pactivity),
            "prior_study_count": count,
            "prior_support_class": "multiple_prior_studies" if count >= 2 else "one_prior_study",
            "heldout_study_bridge_count": int(bridge_counts.get(heldout.document_id, 0)),
            "study_bridge_class": "multi_bridge_study" if int(bridge_counts.get(heldout.document_id, 0)) >= 2 else "single_bridge_study",
            "prior_document_ids": prior_documents,
            "heldout_document_ids": str(heldout.document_id),
            "prior_assay_ids": prior_assays,
            "heldout_assay_ids": heldout.assay_ids,
            "prior_activity_ids": prior_activities,
            "heldout_activity_ids": heldout.activity_ids,
            "component_id": heldout.component_id,
            "is_verified_natural_product": bool(heldout.is_verified_natural_product),
            "history_mode": history_mode,
            "molecular_weight_scaled": float(heldout.molecular_weight_scaled),
            "logp_scaled": float(heldout.logp_scaled),
            **feature_values,
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["heldout_year", "heldout_study_key", "molecule_id"],
    ).reset_index(drop=True) if rows else pd.DataFrame()


class EvidenceConditionedM1:
    """Evidence-conditioned hierarchical model for blind transfer to a new study.

    location_i = historical_median_i + alpha + x_i beta
    log(scale_i) = log_sigma_0 + r_i gamma
    u_study ~ Normal(0, tau^2), residual ~ StudentT(df, 0, scale_i)
    """

    def __init__(
        self, config: M1Config | None = None, *, mean_features: Sequence[str] = MEAN_FEATURES,
        scale_features: Sequence[str] = SCALE_FEATURES, include_study_variance: bool = True,
    ):
        self.config = config or M1Config()
        self.mean_features = tuple(mean_features)
        self.scale_features = tuple(scale_features)
        self.include_study_variance = include_study_variance
        self.feature_means_: dict[str, float] = {}
        self.feature_scales_: dict[str, float] = {}
        self.parameters_: np.ndarray | None = None
        self.optimization_: dict[str, Any] = {}
        nodes, weights = np.polynomial.hermite.hermgauss(self.config.quadrature_nodes)
        self._nodes = nodes
        self._log_weights = np.log(weights) - 0.5 * math.log(math.pi)

    def _fit_transform(self, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        features = tuple(dict.fromkeys((*self.mean_features, *self.scale_features)))
        for feature in features:
            values = frame[feature].to_numpy(float)
            self.feature_means_[feature] = float(np.mean(values))
            scale = float(np.std(values))
            self.feature_scales_[feature] = scale if scale > 1e-8 else 1.0
        return self._transform(frame)

    def _transform(self, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        def matrix(names: Sequence[str]) -> np.ndarray:
            if not names:
                return np.empty((len(frame), 0))
            return np.column_stack([
                (frame[name].to_numpy(float) - self.feature_means_[name]) / self.feature_scales_[name]
                for name in names
            ])
        return matrix(self.mean_features), matrix(self.scale_features)

    def _unpack(self, values: np.ndarray) -> tuple[float, np.ndarray, float, np.ndarray, float]:
        p, q = len(self.mean_features), len(self.scale_features)
        alpha = float(values[0])
        beta = values[1:1 + p]
        log_sigma = float(values[1 + p])
        gamma = values[2 + p:2 + p + q]
        tau = self.config.minimum_scale + math.exp(float(values[-1])) if self.include_study_variance else 0.0
        return alpha, beta, log_sigma, gamma, tau

    def _objective(
        self, values: np.ndarray, frame: pd.DataFrame, mean_x: np.ndarray,
        scale_x: np.ndarray, study_groups: list[np.ndarray],
    ) -> float:
        alpha, beta, log_sigma, gamma, tau = self._unpack(values)
        location = frame.prior_median_pactivity.to_numpy(float) + alpha + mean_x @ beta
        log_scales = np.clip(log_sigma + scale_x @ gamma, -6.0, 3.0)
        scales = self.config.minimum_scale + np.exp(log_scales)
        observed = frame.heldout_median_pactivity.to_numpy(float)
        losses = []
        for indices in study_groups:
            if tau == 0:
                standardized = (observed[indices] - location[indices]) / scales[indices]
                log_likelihood = np.sum(student_t.logpdf(
                    standardized, df=self.config.degrees_of_freedom,
                ) - np.log(scales[indices]))
            else:
                shifts = math.sqrt(2.0) * tau * self._nodes
                standardized = (
                    observed[indices, None] - location[indices, None] - shifts[None, :]
                ) / scales[indices, None]
                node_likelihood = np.sum(
                    student_t.logpdf(standardized, df=self.config.degrees_of_freedom)
                    - np.log(scales[indices, None]), axis=0,
                )
                log_likelihood = logsumexp(self._log_weights + node_likelihood)
            losses.append(-float(log_likelihood) / len(indices))
        penalty = 0.5 * float(np.sum((beta / self.config.mean_coefficient_scale) ** 2))
        penalty += 0.5 * float(np.sum((gamma / self.config.log_scale_coefficient_scale) ** 2))
        return float(np.mean(losses) + penalty / max(len(frame), 1))

    def fit(self, frame: pd.DataFrame) -> "EvidenceConditionedM1":
        required = {
            "prior_median_pactivity", "heldout_median_pactivity", "heldout_study_key",
            *self.mean_features, *self.scale_features,
        }
        if missing := sorted(required - set(frame.columns)):
            raise ValueError(f"Missing M1 columns: {missing}")
        clean = frame.reset_index(drop=True)
        mean_x, scale_x = self._fit_transform(clean)
        residual = clean.heldout_median_pactivity.to_numpy(float) - clean.prior_median_pactivity.to_numpy(float)
        alpha = float(np.median(residual))
        sigma = max(float(np.median(np.abs(residual - alpha)) * 1.4826), 0.1)
        initial = [alpha, *np.zeros(len(self.mean_features)), math.log(sigma), *np.zeros(len(self.scale_features))]
        if self.include_study_variance:
            initial.append(math.log(max(sigma / 2.0, 0.05)))
        groups = [group.index.to_numpy() for _, group in clean.groupby("heldout_study_key", sort=True)]
        bounds = [(-3.0, 3.0)] * len(initial)
        bounds[1 + len(self.mean_features)] = (-6.0, 2.0)
        if self.include_study_variance:
            bounds[-1] = (-7.0, 2.0)
        result = minimize(
            self._objective, np.asarray(initial), args=(clean, mean_x, scale_x, groups),
            method="L-BFGS-B", bounds=bounds,
            options={
                "maxiter": self.config.optimizer_maxiter,
                "maxfun": self.config.optimizer_maxfun,
                "ftol": 1e-10,
                "gtol": 1e-6,
            },
        )
        if not np.isfinite(result.fun):
            raise RuntimeError(f"M1 optimization failed: {result.message}")
        self.parameters_ = result.x
        self.optimization_ = {
            "success": bool(result.success), "message": str(result.message),
            "iterations": int(result.nit), "objective": float(result.fun),
        }
        return self

    def _components(self, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, float]:
        if self.parameters_ is None:
            raise RuntimeError("M1 is not fitted")
        mean_x, scale_x = self._transform(frame)
        alpha, beta, log_sigma, gamma, tau = self._unpack(self.parameters_)
        location = frame.prior_median_pactivity.to_numpy(float) + alpha + mean_x @ beta
        scale = self.config.minimum_scale + np.exp(np.clip(log_sigma + scale_x @ gamma, -6.0, 3.0))
        return location, scale, tau

    def _cdf(self, value: float, location: float, scale: float, tau: float) -> float:
        if tau == 0:
            return float(student_t.cdf((value - location) / scale, df=self.config.degrees_of_freedom))
        shifts = math.sqrt(2.0) * tau * self._nodes
        return float(np.sum(np.exp(self._log_weights) * student_t.cdf(
            (value - location - shifts) / scale, df=self.config.degrees_of_freedom,
        )))

    def predict_distribution(self, frame: pd.DataFrame, levels: Sequence[float] = (0.8, 0.9)) -> dict[str, np.ndarray]:
        location, scale, tau = self._components(frame)
        output: dict[str, np.ndarray] = {
            "prediction": location, "residual_scale": scale,
            "between_study_sd": np.full(len(frame), tau),
            "residual_variance": scale ** 2 * self.config.degrees_of_freedom / (self.config.degrees_of_freedom - 2.0),
            "between_study_variance": np.full(len(frame), tau ** 2),
        }
        if "heldout_median_pactivity" in frame:
            observed = frame.heldout_median_pactivity.to_numpy(float)
            if tau == 0:
                density = student_t.pdf((observed - location) / scale, df=self.config.degrees_of_freedom) / scale
            else:
                shifts = math.sqrt(2.0) * tau * self._nodes
                standardized = (observed[:, None] - location[:, None] - shifts[None, :]) / scale[:, None]
                density = np.exp(logsumexp(
                    self._log_weights[None, :] + student_t.logpdf(
                        standardized, df=self.config.degrees_of_freedom,
                    ) - np.log(scale[:, None]), axis=1,
                ))
            output["nll"] = -np.log(np.maximum(density, 1e-300))
        for level in levels:
            suffix = str(int(level * 100))
            lower, upper = [], []
            tail = (1.0 - level) / 2.0
            for center, row_scale in zip(location, scale):
                variance = row_scale ** 2 * self.config.degrees_of_freedom / (self.config.degrees_of_freedom - 2.0) + tau ** 2
                radius = max(12.0 * math.sqrt(variance), 5.0)
                quantiles = [brentq(
                    lambda value: self._cdf(value, float(center), float(row_scale), tau) - probability,
                    float(center - radius), float(center + radius),
                ) for probability in (tail, 1.0 - tail)]
                lower.append(quantiles[0]); upper.append(quantiles[1])
            output[f"lower_{suffix}"] = np.asarray(lower)
            output[f"upper_{suffix}"] = np.asarray(upper)
        return output

    def parameters(self) -> dict[str, Any]:
        if self.parameters_ is None:
            raise RuntimeError("M1 is not fitted")
        alpha, beta, log_sigma, gamma, tau = self._unpack(self.parameters_)
        return {
            "alpha": alpha,
            "mean_coefficients": dict(zip(self.mean_features, map(float, beta))),
            "base_residual_scale": float(self.config.minimum_scale + math.exp(log_sigma)),
            "log_scale_coefficients": dict(zip(self.scale_features, map(float, gamma))),
            "between_study_sd": float(tau),
            "feature_means": self.feature_means_, "feature_scales": self.feature_scales_,
            "degrees_of_freedom": self.config.degrees_of_freedom,
            "include_study_variance": self.include_study_variance,
            "optimization": self.optimization_,
        }


def _t_distribution(prediction: np.ndarray, outcome: np.ndarray, scale: float, df: float) -> dict[str, np.ndarray]:
    result = {
        "prediction": prediction,
        "nll": -student_t.logpdf((outcome - prediction) / scale, df=df) + math.log(scale),
        "residual_scale": np.full(len(prediction), scale),
        "between_study_sd": np.zeros(len(prediction)),
        "residual_variance": np.full(len(prediction), scale ** 2 * df / (df - 2.0)),
        "between_study_variance": np.zeros(len(prediction)),
    }
    for level in (0.8, 0.9):
        critical = student_t.ppf((1.0 + level) / 2.0, df=df)
        suffix = str(int(level * 100))
        result[f"lower_{suffix}"] = prediction - critical * scale
        result[f"upper_{suffix}"] = prediction + critical * scale
    return result


def _xgb_model(seed: int) -> XGBRegressor:
    return XGBRegressor(
        objective="reg:squarederror", n_estimators=120, max_depth=2,
        learning_rate=0.03, min_child_weight=5, subsample=0.8,
        colsample_bytree=0.8, reg_lambda=10.0, reg_alpha=0.1,
        random_state=seed, n_jobs=1, verbosity=0,
    )


def _xgb_distribution(train: pd.DataFrame, test: pd.DataFrame, config: M1Config, seed: int) -> dict[str, np.ndarray]:
    x_train = train.loc[:, ALL_TABULAR_FEATURES].to_numpy(float)
    delta = train.heldout_median_pactivity.to_numpy(float) - train.prior_median_pactivity.to_numpy(float)
    groups = train.heldout_study_key.astype(str).to_numpy()
    unique_groups = np.unique(groups)
    cross_fitted = np.full(len(train), np.nan)
    if len(unique_groups) >= 3:
        splitter = GroupKFold(n_splits=min(5, len(unique_groups)))
        for fold, (fit_idx, valid_idx) in enumerate(splitter.split(x_train, delta, groups)):
            model = _xgb_model(seed + fold).fit(x_train[fit_idx], delta[fit_idx])
            cross_fitted[valid_idx] = model.predict(x_train[valid_idx])
    if np.isnan(cross_fitted).any():
        cross_fitted[np.isnan(cross_fitted)] = float(np.median(delta))
    scale = _student_scale(delta - cross_fitted, config.degrees_of_freedom)
    model = _xgb_model(seed).fit(x_train, delta)
    prediction = test.prior_median_pactivity.to_numpy(float) + model.predict(
        test.loc[:, ALL_TABULAR_FEATURES].to_numpy(float),
    )
    return _t_distribution(
        prediction, test.heldout_median_pactivity.to_numpy(float), scale, config.degrees_of_freedom,
    )


def _identity(test: pd.DataFrame, scheme: str, fold: int) -> pd.DataFrame:
    columns = [
        "target_chembl_id", "molecule_id", "heldout_study_key", "heldout_year",
        "prior_median_pactivity", "heldout_median_pactivity", "prior_study_count",
        "prior_support_class", "heldout_study_bridge_count", "study_bridge_class",
        "prior_document_ids", "heldout_document_ids", "prior_assay_ids", "heldout_assay_ids",
        "prior_activity_ids", "heldout_activity_ids", "component_id",
        "is_verified_natural_product", "history_mode", *ALL_TABULAR_FEATURES,
    ]
    output = test.loc[:, list(dict.fromkeys(columns))].copy().reset_index(drop=True)
    output.insert(0, "prediction_id", [f"{scheme}:{fold}:{index}" for index in range(len(output))])
    return output


def _fit_fold(
    train_edges: pd.DataFrame, test_edges: pd.DataFrame, config: M1Config, *,
    scheme: str, fold: int, history_mode: str, temporal: bool, seed: int,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    train = construct_feature_events(train_edges, history_mode="leave_publication_out_development")
    test = construct_feature_events(train_edges, test_edges, history_mode=history_mode)
    if len(train) < 8 or train.heldout_study_key.nunique() < 4 or test.empty:
        return pd.DataFrame(), []
    assert_fold_provenance(train_edges, test, temporal=temporal)
    m0_config = M0Config(
        degrees_of_freedom=config.degrees_of_freedom,
        quadrature_nodes=config.quadrature_nodes,
        minimum_scale=config.minimum_scale,
        optimizer_maxiter=config.optimizer_maxiter,
    )
    student = HierarchicalStudentTM0(m0_config, include_study_variance=False).fit(train)
    m0 = HierarchicalStudentTM0(m0_config, include_study_variance=True).fit(train)
    no_protocol_scale = tuple(feature for feature in SCALE_FEATURES if feature not in PROTOCOL_SCALE_FEATURES)
    models = {
        "m1_no_protocol": EvidenceConditionedM1(config, scale_features=no_protocol_scale).fit(train),
        "m1_homoscedastic": EvidenceConditionedM1(config, scale_features=()).fit(train),
        "m1_no_study_variance": EvidenceConditionedM1(config, include_study_variance=False).fit(train),
        "phytotrust_m1": EvidenceConditionedM1(config).fit(train),
    }
    outcome = test.heldout_median_pactivity.to_numpy(float)
    median_prediction = test.prior_median_pactivity.to_numpy(float)
    median_scale = _student_scale(
        train.heldout_median_pactivity.to_numpy(float) - train.prior_median_pactivity.to_numpy(float),
        config.degrees_of_freedom,
    )
    predictions = {
        "median_potency": _t_distribution(median_prediction, outcome, median_scale, config.degrees_of_freedom),
        "student_t_regression": student.predict_distribution(test, levels=(0.8, 0.9)),
        "xgboost": _xgb_distribution(train, test, config, seed + fold),
        "phytotrust_m0": m0.predict_distribution(test, levels=(0.8, 0.9)),
        **{name: model.predict_distribution(test) for name, model in models.items()},
    }
    identity = _identity(test, scheme, fold)
    rows = []
    for method, values in predictions.items():
        result = identity.copy()
        result.insert(0, "method", method); result.insert(0, "fold", fold)
        result.insert(0, "validation_scheme", scheme)
        for name, value in values.items():
            result[name] = value
        rows.append(result)
    parameters = []
    for name, model in models.items():
        values = model.parameters()
        parameters.append({
            "validation_scheme": scheme, "fold": fold, "method": name,
            "training_events": len(train), "training_studies": train.heldout_study_key.nunique(),
            "test_events": len(test), "test_studies": test.heldout_study_key.nunique(),
            "alpha": values["alpha"], "base_residual_scale": values["base_residual_scale"],
            "between_study_sd": values["between_study_sd"],
            "mean_coefficients": json.dumps(values["mean_coefficients"], sort_keys=True),
            "log_scale_coefficients": json.dumps(values["log_scale_coefficients"], sort_keys=True),
            "optimizer_success": values["optimization"]["success"],
            "optimizer_message": values["optimization"]["message"],
        })
    return pd.concat(rows, ignore_index=True), parameters


def cross_validate(edges: pd.DataFrame, config: M1Config, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    results: list[pd.DataFrame] = []
    parameters: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    repeated = edges.groupby("molecule_id").document_id.transform("nunique") >= 2
    for fold, study in enumerate(sorted(edges.loc[repeated, "document_id"].unique())):
        result, params = _fit_fold(
            edges[edges.document_id != study], edges[edges.document_id == study], config,
            scheme="study_grouped", fold=fold,
            history_mode="other_development_studies_not_time_restricted",
            temporal=False, seed=seed,
        )
        if result.empty:
            skipped.append({"scheme": "study_grouped", "group": study, "reason": "insufficient_connected_training_data"})
        else:
            results.append(result); parameters.extend(params)
    temporal_fold = 0
    for year in sorted(edges.document_year.unique()):
        train_edges = edges[edges.document_year < year]
        if train_edges.empty:
            continue
        result, params = _fit_fold(
            train_edges, edges[edges.document_year == year], config,
            scheme="forward_temporal", fold=temporal_fold,
            history_mode="strictly_earlier_publications", temporal=True, seed=seed,
        )
        if result.empty:
            skipped.append({"scheme": "forward_temporal", "group": int(year), "reason": "no_eligible_transfer_events_or_connected_training_data"})
        else:
            results.append(result); parameters.extend(params)
        temporal_fold += 1
    if not results:
        raise RuntimeError("No M1 validation folds were estimable")
    return pd.concat(results, ignore_index=True), pd.DataFrame(parameters), skipped


def _metric(group: pd.DataFrame) -> dict[str, Any]:
    observed = group.heldout_median_pactivity.to_numpy(float)
    prediction = group.prediction.to_numpy(float)
    result = {
        "n": int(len(group)), "compounds": int(group.molecule_id.nunique()),
        "studies": int(group.heldout_study_key.nunique()),
        "mae": float(mean_absolute_error(observed, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(observed, prediction))),
        "nll": float(group.nll.mean()),
    }
    calibration = []
    for level in (80, 90):
        coverage = float(((observed >= group[f"lower_{level}"]) & (observed <= group[f"upper_{level}"])).mean())
        result[f"coverage_{level}"] = coverage
        result[f"mean_width_{level}"] = float((group[f"upper_{level}"] - group[f"lower_{level}"]).mean())
        calibration.append(abs(coverage - level / 100.0))
    result["mean_absolute_calibration_error"] = float(np.mean(calibration))
    return result


def evaluate(predictions: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    scopes = {
        "all": lambda value: value,
        "natural_products": lambda value: value[value.is_verified_natural_product],
        "one_prior_study": lambda value: value[value.prior_support_class == "one_prior_study"],
        "multiple_prior_studies": lambda value: value[value.prior_support_class == "multiple_prior_studies"],
    }
    rows, nested = [], {}
    for (scheme, method), group in predictions.groupby(["validation_scheme", "method"], sort=True):
        for scope, selector in scopes.items():
            selected = selector(group)
            values = _metric(selected) if len(selected) else {}
            rows.append({"validation_scheme": scheme, "scope": scope, "method": method, **values})
            nested.setdefault(scheme, {}).setdefault(scope, {})[method] = values
    return pd.DataFrame(rows), nested


def cluster_bootstrap(
    predictions: pd.DataFrame, *, iterations: int, seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    records, nested = [], {}
    for scheme in sorted(predictions.validation_scheme.unique()):
        scheme_frame = predictions[predictions.validation_scheme == scheme]
        for scope, frame in (
            ("all", scheme_frame),
            ("natural_products", scheme_frame[scheme_frame.is_verified_natural_product]),
        ):
            m1 = frame[frame.method == "phytotrust_m1"].set_index("prediction_id")
            studies = sorted(m1.heldout_study_key.unique())
            if len(studies) < 2:
                continue
            nested.setdefault(scheme, {}).setdefault(scope, {})
            clusters = {study: m1.index[m1.heldout_study_key == study].tolist() for study in studies}
            for comparator in METHODS[:-1]:
                other = frame[frame.method == comparator].set_index("prediction_id").loc[m1.index]
                truth = m1.heldout_median_pactivity.to_numpy(float)
                point = {
                    "mae": float(np.mean(np.abs(truth - m1.prediction)) - np.mean(np.abs(truth - other.prediction))),
                    "rmse": float(np.sqrt(np.mean((truth - m1.prediction) ** 2)) - np.sqrt(np.mean((truth - other.prediction) ** 2))),
                    "nll": float(m1.nll.mean() - other.nll.mean()),
                    "calibration": float(_metric(m1.reset_index())["mean_absolute_calibration_error"] - _metric(other.reset_index())["mean_absolute_calibration_error"]),
                }
                samples = {key: [] for key in point}
                for _ in range(iterations):
                    chosen = rng.choice(studies, len(studies), replace=True)
                    identifiers = [item for study in chosen for item in clusters[study]]
                    left, right = m1.loc[identifiers], other.loc[identifiers]
                    y = left.heldout_median_pactivity.to_numpy(float)
                    samples["mae"].append(float(np.mean(np.abs(y - left.prediction)) - np.mean(np.abs(y - right.prediction))))
                    samples["rmse"].append(float(np.sqrt(np.mean((y - left.prediction) ** 2)) - np.sqrt(np.mean((y - right.prediction) ** 2))))
                    samples["nll"].append(float(left.nll.mean() - right.nll.mean()))
                    samples["calibration"].append(float(_metric(left.reset_index())["mean_absolute_calibration_error"] - _metric(right.reset_index())["mean_absolute_calibration_error"]))
                ci = {key: [float(value) for value in np.quantile(sample, [0.025, 0.975])] for key, sample in samples.items()}
                nested[scheme][scope][comparator] = {
                    "point": point, "ci95": ci,
                    "ci_excludes_zero": {key: bool(value[1] < 0 or value[0] > 0) for key, value in ci.items()},
                }
                records.append({
                    "validation_scheme": scheme, "scope": scope,
                    "comparison": f"phytotrust_m1-minus-{comparator}",
                    **{f"delta_{key}": value for key, value in point.items()},
                    **{f"delta_{key}_ci95": json.dumps(value) for key, value in ci.items()},
                })
    return pd.DataFrame(records), nested


def protocol_coverage(edges: pd.DataFrame) -> dict[str, Any]:
    return {
        "compound_study_edges": int(len(edges)),
        "assays": int(edges.protocol_assay_count.sum()),
        "mean_direct_field_completeness": float(edges.protocol_direct_completeness.mean()),
        "mean_protocol_missing_fraction": float(edges.protocol_missing_fraction.mean()),
        "edges_with_any_extracted_protocol_category": int((edges.protocol_category_count > 0).sum()),
        "fraction_with_any_extracted_protocol_category": float((edges.protocol_category_count > 0).mean()),
        "median_protocol_category_count": float(edges.protocol_category_count.median()),
    }


def _render_report(results: dict[str, Any], metrics_table: pd.DataFrame) -> str:
    dataset, coverage = results["dataset"], results["protocol_coverage"]
    lines = [
        "# PhytoTRUST-M1 development-only scientific evaluation", "",
        "## Scope", "",
        f"M1 used only the CHEMBL222 development artifact ({dataset['raw_measurements']} measurements, {dataset['compounds']} compounds, {dataset['studies']} publications) and full ChEMBL metadata for {coverage['assays']} referenced assay instances. The future-evaluation history, index, and quarantined labels were not opened.", "",
        "The complete frozen target snapshot had already been exhausted. No post-outcome compatibility relaxation or outcome-guided measurement acquisition was performed. The targeted data improvement was development-only assay-metadata retrieval and literal protocol extraction.", "",
        "## Method", "",
        "M1 predicts a new-study measurement as the historical compound median plus a regularized evidence-conditioned correction and a publication random intercept. Its Student-t residual scale is conditioned on historical dispersion, support, protocol missingness/diversity, and history span. The unseen publication effect is integrated out. Held-out protocol fields are never predictors.", "",
        "The individual ingredients are established statistical ideas. The tested contribution is their timing-safe integration for new-study evidence transferability with explicit evidence-conditioned uncertainty; this experiment does not establish novelty.", "",
        "## Protocol metadata coverage", "",
        f"Mean direct-field completeness was {coverage['mean_direct_field_completeness']:.3f}; mean missing fraction across the prespecified protocol fields was {coverage['mean_protocol_missing_fraction']:.3f}. Literal protocol categories were found for {coverage['edges_with_any_extracted_protocol_category']}/{coverage['compound_study_edges']} compound-publication edges. Unknown values remained unknown.", "",
        "## Results", "",
        "| Validation | Scope | Method | n | MAE | RMSE | NLL | 80% coverage / width | 90% coverage / width | calibration error |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics_table[metrics_table.scope.isin(["all", "natural_products"])].itertuples():
        if pd.isna(getattr(row, "n", np.nan)):
            continue
        lines.append(
            f"| {row.validation_scheme} | {row.scope} | {row.method} | {int(row.n)} | {row.mae:.3f} | {row.rmse:.3f} | {row.nll:.3f} | {row.coverage_80:.3f} / {row.mean_width_80:.3f} | {row.coverage_90:.3f} / {row.mean_width_90:.3f} | {row.mean_absolute_calibration_error:.3f} |"
        )
    sg = results["metrics"]["study_grouped"]["all"]
    ft = results["metrics"]["forward_temporal"]["all"]
    strongest_sg = min((name for name in METHODS[:-1]), key=lambda name: sg[name]["nll"])
    strongest_ft = min((name for name in METHODS[:-1]), key=lambda name: ft[name]["nll"])
    comparisons = results["study_cluster_bootstrap_differences"]
    lines += [
        "", "## Primary interpretation", "",
        f"Study-grouped M1 NLL was {sg['phytotrust_m1']['nll']:.3f}; the strongest comparator was {strongest_sg} at {sg[strongest_sg]['nll']:.3f}. Forward-temporal M1 NLL was {ft['phytotrust_m1']['nll']:.3f}; the strongest comparator was {strongest_ft} at {ft[strongest_ft]['nll']:.3f}.", "",
        f"M1-minus-{strongest_sg} study-grouped NLL 95% CI: {comparisons['study_grouped']['all'][strongest_sg]['ci95']['nll']}. M1-minus-{strongest_ft} forward-temporal NLL 95% CI: {comparisons['forward_temporal']['all'][strongest_ft]['ci95']['nll']}. Negative values favor M1.", "",
        "## Ablations", "",
        f"Study-grouped NLL: full M1 {sg['phytotrust_m1']['nll']:.3f}; without historical protocol features {sg['m1_no_protocol']['nll']:.3f}; homoscedastic M1 {sg['m1_homoscedastic']['nll']:.3f}; without study variance {sg['m1_no_study_variance']['nll']:.3f}; M0 {sg['phytotrust_m0']['nll']:.3f}.", "",
        f"Forward-temporal NLL: full M1 {ft['phytotrust_m1']['nll']:.3f}; without protocol {ft['m1_no_protocol']['nll']:.3f}; homoscedastic {ft['m1_homoscedastic']['nll']:.3f}; without study variance {ft['m1_no_study_variance']['nll']:.3f}; M0 {ft['phytotrust_m0']['nll']:.3f}.", "",
        "An apparent gain is attributable to the proposed protocol-conditioned scale only if full M1 improves over `m1_no_protocol` with uncertainty excluding zero in both validation schemes. The same logic applies to heteroscedasticity and study variance.", "",
        "## Natural products and connectivity", "",
        f"Only {dataset['verified_natural_product_transfer_events']} connected natural-product transfer events were available in the full development graph. NP results are exploratory. The multiple-prior-study stratum contains only {results['connectivity']['study_grouped_multiple_prior_events']} study-grouped and {results['connectivity']['forward_temporal_multiple_prior_events']} temporal events, so improved performance cannot be attributed reliably to deeper compound connectivity.", "",
        "## Limitations", "",
        "- M1 predicts measurement transferability, not true bioactivity, efficacy, safety, or clinical utility.",
        "- Publication is an imperfect proxy for independent study and laboratory.",
        "- ChEMBL protocol descriptions are retrospective records. Only historical descriptions were used, but their original pre-experiment availability is not proven.",
        "- Protocol fields are sparse and literal text extraction cannot recover unreported conditions.",
        "- Leave-publication-out compound residuals are crossed and dependent; publication bootstrap intervals remain approximate.",
        "- The study-grouped analysis is retrospective; only forward-temporal validation obeys a strict time filtration.",
        "- One target, 29 repeated compounds, and very small NP strata prevent broad generalization.",
        "", "## Publication-readiness verdict", "",
        results["scientific_verdict"],
    ]
    return "\n".join(lines) + "\n"


def _plot_metric(metrics: pd.DataFrame, path: Path, metric: str, title: str) -> None:
    rows = metrics[(metrics.scope == "all") & metrics.method.isin(METHODS)]
    colors = ["#0ea5e9", "#a78bfa", "#f59e0b", "#22c55e", "#64748b", "#14b8a6", "#ec4899", "#f8fafc"]
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="720" viewBox="0 0 1280 720">',
        '<rect width="1280" height="720" fill="#07111f"/>',
        f'<text x="55" y="45" fill="#f8fafc" font-family="sans-serif" font-size="24">{title}</text>',
        '<text x="55" y="70" fill="#94a3b8" font-family="sans-serif" font-size="13">Development-only; lower is better.</text>',
    ]
    for panel, scheme in enumerate(("study_grouped", "forward_temporal")):
        subset = rows[rows.validation_scheme == scheme]
        x0 = 55 + panel * 610
        maximum = max(float(subset[metric].max()), 0.001)
        pieces.append(f'<text x="{x0}" y="115" fill="#e2e8f0" font-family="sans-serif" font-size="18">{scheme.replace("_", " ").title()}</text>')
        for index, (method, color) in enumerate(zip(METHODS, colors)):
            match = subset[subset.method == method]
            if match.empty:
                continue
            value = float(match.iloc[0][metric]); y = 145 + index * 62
            width = 330 * value / maximum
            pieces.append(f'<text x="{x0}" y="{y+16}" fill="#cbd5e1" font-family="sans-serif" font-size="12">{method}</text>')
            pieces.append(f'<rect x="{x0+190}" y="{y}" width="{width:.1f}" height="24" rx="4" fill="{color}"/>')
            pieces.append(f'<text x="{x0+197+width:.1f}" y="{y+17}" fill="#f8fafc" font-family="sans-serif" font-size="12">{value:.3f}</text>')
    pieces.append("</svg>")
    path.write_text("\n".join(pieces), encoding="utf-8")


def run_experiment(
    dataset: Path, protocol_features_path: Path, output_dir: Path, *,
    bootstrap_iterations: int = 1000, seed: int = SEED, config: M1Config | None = None,
) -> dict[str, Any]:
    config = config or M1Config()
    measurements, dataset_manifest = load_connected_development(dataset)
    protocol, protocol_manifest = load_protocol_features(
        protocol_features_path, set(measurements.assay_id.astype(str)),
    )
    edges = enrich_edges(measurements, protocol)
    predictions, fold_parameters, skipped = cross_validate(edges, config, seed)
    metrics_table, metrics = evaluate(predictions)
    bootstrap_table, bootstrap = cluster_bootstrap(
        predictions, iterations=bootstrap_iterations, seed=seed,
    )
    final_events = construct_feature_events(edges, history_mode="leave_publication_out_development")
    final_model = EvidenceConditionedM1(config).fit(final_events)
    coverage = protocol_coverage(edges)
    dataset_manifest.update({
        "compound_study_edges": int(len(edges)),
        "connected_transfer_events": int(len(final_events)),
        "repeated_compounds": int(final_events.molecule_id.nunique()),
        "connected_studies": int(final_events.heldout_study_key.nunique()),
        "verified_natural_product_transfer_events": int(final_events.is_verified_natural_product.sum()),
        "protocol_artifact": protocol_manifest,
    })
    sg, ft = metrics["study_grouped"]["all"], metrics["forward_temporal"]["all"]
    strongest_sg = min(METHODS[:-1], key=lambda method: sg[method]["nll"])
    strongest_ft = min(METHODS[:-1], key=lambda method: ft[method]["nll"])
    sg_ci = bootstrap["study_grouped"]["all"][strongest_sg]["ci95"]["nll"]
    ft_ci = bootstrap["forward_temporal"]["all"][strongest_ft]["ci95"]["nll"]
    no_np_degradation = (
        metrics["study_grouped"]["natural_products"]["phytotrust_m1"]["nll"]
        <= metrics["study_grouped"]["natural_products"][strongest_sg]["nll"]
        and metrics["forward_temporal"]["natural_products"]["phytotrust_m1"]["nll"]
        <= metrics["forward_temporal"]["natural_products"][strongest_ft]["nll"]
    )
    supported = bool(
        sg_ci[1] < 0 and ft_ci[1] < 0
        and sg["phytotrust_m1"]["mean_absolute_calibration_error"] <= sg[strongest_sg]["mean_absolute_calibration_error"]
        and ft["phytotrust_m1"]["mean_absolute_calibration_error"] <= ft[strongest_ft]["mean_absolute_calibration_error"]
        and no_np_degradation
    )
    verdict = (
        "M1 met the frozen development-only support rule against the strongest baselines. It is eligible for one genuinely untouched external evaluation, but this is not proof of novelty or journal-level validity."
        if supported else
        "M1 did not meet the frozen development-only support rule. The evidence-conditioned formulation remains a valid negative feasibility result, but the available data do not support a claim of superiority, methodological novelty, or IEEE journal readiness. M1 must not be promoted to the quarantined final test."
    )
    results = {
        "generated_at": datetime.now(UTC).isoformat(), "model_version": MODEL_VERSION,
        "seed": seed, "configuration": asdict(config), "dataset": dataset_manifest,
        "protocol_coverage": coverage,
        "validation": {
            "study_grouped": "retrospective leave-one-publication-out",
            "forward_temporal": "strictly earlier publication years only",
            "heldout_protocol_used_as_predictor": False,
            "heldout_publication_assay_activity_outcome_excluded": True,
            "quarantined_artifacts_loaded": False,
            "skipped_folds": skipped,
        },
        "metrics": metrics,
        "study_cluster_bootstrap_differences": bootstrap,
        "connectivity": {
            "study_grouped_multiple_prior_events": sg["phytotrust_m1"]["n"] - metrics["study_grouped"]["one_prior_study"]["phytotrust_m1"]["n"],
            "forward_temporal_multiple_prior_events": ft["phytotrust_m1"]["n"] - metrics["forward_temporal"]["one_prior_study"]["phytotrust_m1"]["n"],
        },
        "final_parameters": final_model.parameters(),
        "strongest_baselines": {"study_grouped": strongest_sg, "forward_temporal": strongest_ft},
        "support_rule_passed": supported,
        "scientific_verdict": verdict,
        "scientific_status": "development_only_m1_not_independently_validated",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "predictions.csv", index=False)
    metrics_table.to_csv(output_dir / "metrics.csv", index=False)
    fold_parameters.to_csv(output_dir / "fold_parameters.csv", index=False)
    bootstrap_table.to_csv(output_dir / "bootstrap_differences.csv", index=False)
    edges.to_csv(output_dir / "enriched_compound_study_edges.csv", index=False)
    final_events.to_csv(output_dir / "m1_connected_events.csv", index=False)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir / "model_config.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "configuration": asdict(config),
        "mean_features": MEAN_FEATURES, "scale_features": SCALE_FEATURES,
        "equations": {
            "location": "mu_i = historical_median_i + alpha + beta^T x_i",
            "scale": "log sigma_i = log sigma_0 + gamma^T r_i",
            "measurement": "Y_is = mu_i + u_s + StudentT(df, 0, sigma_i)",
            "study_effect": "u_s ~ Normal(0, tau^2); integrate u_new for blind transfer",
        },
    }, indent=2), encoding="utf-8")
    (output_dir / "model_parameters.json").write_text(json.dumps(final_model.parameters(), indent=2), encoding="utf-8")
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_render_report(results, metrics_table), encoding="utf-8")
    _plot_metric(metrics_table, output_dir / "nll_comparison.svg", "nll", "PhytoTRUST-M1 probabilistic comparison")
    _plot_metric(metrics_table, output_dir / "calibration_comparison.svg", "mean_absolute_calibration_error", "PhytoTRUST-M1 calibration comparison")
    (output_dir / "command.txt").write_text(
        "python -m app.phytotrust_m1 run --dataset ../experiments/phytotrust_connected_dataset/development_measurements.csv --protocol-features ../data/phytotrust_m1/CHEMBL222.assay_protocol_features.csv --output-dir ../experiments/phytotrust_m1 --bootstrap-iterations 1000 --seed 2026\n",
        encoding="utf-8",
    )
    artifact_names = sorted(path.name for path in output_dir.iterdir() if path.is_file() and path.name != "manifest.json")
    (output_dir / "manifest.json").write_text(json.dumps({
        "model_version": MODEL_VERSION, "development_source_sha256": dataset_manifest["source_sha256"],
        "protocol_source_sha256": protocol_manifest["sha256"],
        "files": {name: _sha256(output_dir / name) for name in artifact_names},
        "quarantined_artifacts_loaded": False,
    }, indent=2), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="PhytoTRUST-M1 evidence-conditioned transferability experiment")
    subparsers = parser.add_subparsers(dest="command", required=True)
    fetch = subparsers.add_parser("fetch-assays")
    fetch.add_argument("--dataset", type=Path, required=True)
    fetch.add_argument("--raw-output", type=Path, required=True)
    fetch.add_argument("--feature-output", type=Path, required=True)
    fetch.add_argument("--batch-size", type=int, default=50)
    fetch.add_argument("--max-retries", type=int, default=5)
    run = subparsers.add_parser("run")
    run.add_argument("--dataset", type=Path, required=True)
    run.add_argument("--protocol-features", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--bootstrap-iterations", type=int, default=1000)
    run.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    if args.command == "fetch-assays":
        result = asyncio.run(fetch_development_assays(
            args.dataset, args.raw_output, args.feature_output,
            batch_size=args.batch_size, max_retries=args.max_retries,
        ))
    else:
        result = run_experiment(
            args.dataset, args.protocol_features, args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations, seed=args.seed,
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
