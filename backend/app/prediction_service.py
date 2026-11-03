from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .features import molecular_descriptors
from .ml import predict_activity


class ModelArtifactError(RuntimeError):
    """The configured inference artifact is absent, corrupt, or incompatible."""


class PredictionError(RuntimeError):
    """The loaded model could not produce a prediction."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def default_model_directory() -> Path:
    return Path(__file__).resolve().parents[2] / "experiments" / "ml_activity_chembl335"


def configured_model_directory() -> Path:
    configured = os.getenv("PHYTOTRUTH_MODEL_DIR")
    return Path(configured).expanduser().resolve() if configured else default_model_directory()


@dataclass(frozen=True)
class ActivityPredictor:
    model: Any
    config: dict[str, Any]
    reference_fingerprints: np.ndarray

    def predict(self, smiles: str) -> dict[str, Any]:
        bundle = {
            "model": self.model,
            "feature_config": self.config["feature_config"],
            "uncertainty": self.config["uncertainty"],
            "applicability_domain": {
                "similarity_threshold": self.config["applicability_domain"]["similarity_threshold"],
                "reference_packed_fingerprints": self.reference_fingerprints,
            },
            "training": self.config["training"],
        }
        try:
            result = predict_activity(bundle, smiles)
            descriptors = molecular_descriptors(result["canonical_smiles"])
        except ValueError:
            raise
        except Exception as exc:  # model/runtime failures are distinct from invalid input
            raise PredictionError("The activity model failed during inference") from exc

        inside = bool(result["inside_applicability_domain"])
        limitations = list(self.config["prediction_semantics"]["limitations"])
        if not inside:
            limitations.insert(0, "OUTSIDE APPLICABILITY DOMAIN: do not treat this prediction as reliable.")
        return {
            "canonical_smiles": result["canonical_smiles"],
            "predicted_pchembl": result["predicted_pactivity"],
            "prediction_interval": {
                "lower": result["prediction_interval_90"][0],
                "upper": result["prediction_interval_90"][1],
                "confidence_level": float(self.config["prediction_semantics"]["confidence_level"]),
                "method": result["interval_method"],
            },
            "applicability_domain": {
                "inside": inside,
                "status": "IN_DOMAIN_RESEARCH_ONLY" if inside else "OUTSIDE_DOMAIN_DO_NOT_RELY",
                "nearest_training_tanimoto": result["nearest_training_tanimoto"],
                "threshold": result["applicability_domain_threshold"],
                "warning": result["warning"],
            },
            "descriptors": descriptors,
            "model": {
                "name": self.config["model"]["name"],
                "version": self.config["model"]["version"],
                "format": self.config["model"]["format"],
                "target": self.config["training"]["target"],
                "training_structures": self.config["training"]["n_structures"],
                "training_scaffolds": self.config["training"]["n_scaffolds"],
            },
            "kind": "predicted",
            "limitations": limitations,
        }


@lru_cache(maxsize=4)
def load_predictor(model_directory: str) -> ActivityPredictor:
    """Load and validate native model artifacts once per resolved directory."""
    try:
        import xgboost
    except ImportError as exc:  # pragma: no cover - installation failure
        raise ModelArtifactError("XGBoost runtime is unavailable") from exc

    root = Path(model_directory)
    config_path = root / "inference_config.json"
    if not config_path.is_file():
        raise ModelArtifactError(f"Missing model configuration: {config_path}")
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        model_path = root / config["model"]["file"]
        fingerprints_path = root / config["applicability_domain"]["fingerprints_file"]
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ModelArtifactError("Inference configuration is invalid") from exc
    if config.get("artifact_version") != 2 or config["model"].get("format") != "xgboost-json":
        raise ModelArtifactError("Unsupported inference artifact format")
    if not model_path.is_file() or not fingerprints_path.is_file():
        raise ModelArtifactError("Native XGBoost model or applicability fingerprints are missing")
    if _sha256(model_path) != config["model"].get("sha256"):
        raise ModelArtifactError("XGBoost model checksum mismatch")
    if _sha256(fingerprints_path) != config["applicability_domain"].get("fingerprints_sha256"):
        raise ModelArtifactError("Applicability fingerprint checksum mismatch")

    try:
        model = xgboost.XGBRegressor()
        model.load_model(model_path)
        # Bound OpenMP work per request so concurrent API calls do not oversubscribe
        # every CPU core. This does not alter the fitted trees or predictions.
        model.set_params(n_jobs=1)
        model.get_booster().set_param({"nthread": 1})
        with np.load(fingerprints_path, allow_pickle=False) as archive:
            fingerprints = np.asarray(archive["reference_packed_fingerprints"], dtype=np.uint8)
    except Exception as exc:
        raise ModelArtifactError("Native inference artifacts could not be loaded") from exc
    feature_count = int(config["feature_config"]["feature_count"])
    n_bits = int(config["feature_config"]["n_bits"])
    expected_width = (n_bits + 7) // 8
    if int(model.n_features_in_) != feature_count:
        raise ModelArtifactError("Model feature count does not match feature configuration")
    if fingerprints.ndim != 2 or fingerprints.shape != (
        int(config["applicability_domain"]["reference_count"]), expected_width
    ):
        raise ModelArtifactError("Applicability fingerprints do not match their configuration")
    fingerprints.setflags(write=False)
    return ActivityPredictor(model=model, config=config, reference_fingerprints=fingerprints)


def get_predictor() -> ActivityPredictor:
    return load_predictor(str(configured_model_directory()))


def clear_model_cache() -> None:
    load_predictor.cache_clear()
