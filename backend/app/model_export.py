from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_inference_artifacts(source: Path, output_dir: Path) -> dict[str, Any]:
    """Convert the trained joblib bundle to supported XGBoost JSON + safe data files.

    This is a format conversion only. It never fits or mutates model parameters.
    """
    import joblib
    import xgboost

    if not source.is_file():
        raise FileNotFoundError(source)
    bundle = joblib.load(source)
    required = {"model", "feature_config", "uncertainty", "applicability_domain", "training"}
    if not isinstance(bundle, dict) or not required.issubset(bundle):
        raise ValueError("Source is not a valid PhytoTruth model bundle")
    if bundle["training"].get("model_name") != "xgboost":
        raise ValueError("Native export currently requires the selected XGBoost model")

    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / "model.xgb.json"
    fingerprints_path = output_dir / "reference_fingerprints.npz"
    config_path = output_dir / "inference_config.json"
    bundle["model"].save_model(model_path)
    fingerprints = np.asarray(
        bundle["applicability_domain"]["reference_packed_fingerprints"],
        dtype=np.uint8,
    )
    if fingerprints.ndim != 2 or not len(fingerprints):
        raise ValueError("Model bundle has invalid reference fingerprints")
    np.savez_compressed(fingerprints_path, reference_packed_fingerprints=fingerprints)

    model_hash = _sha256(model_path)
    fingerprint_hash = _sha256(fingerprints_path)
    version = f"chembl335-xgboost-{model_hash[:12]}"
    config = {
        "artifact_version": 2,
        "exported_at": datetime.now(UTC).isoformat(),
        "model": {
            "name": "xgboost",
            "format": "xgboost-json",
            "file": model_path.name,
            "sha256": model_hash,
            "version": version,
            "exported_with_xgboost": xgboost.__version__,
        },
        "feature_config": bundle["feature_config"],
        "uncertainty": bundle["uncertainty"],
        "applicability_domain": {
            "similarity_threshold": bundle["applicability_domain"]["similarity_threshold"],
            "method": bundle["applicability_domain"]["method"],
            "fingerprints_file": fingerprints_path.name,
            "fingerprints_sha256": fingerprint_hash,
            "reference_count": int(fingerprints.shape[0]),
            "packed_bytes_per_fingerprint": int(fingerprints.shape[1]),
        },
        "training": bundle["training"],
        "prediction_semantics": {
            "output": "target-specific experimental pActivity estimate, reported as predicted pChEMBL",
            "confidence_level": bundle["uncertainty"]["coverage_target"],
            "limitations": [
                "Research prioritization only; not confirmed experimental activity.",
                "Not evidence of drug efficacy, clinical benefit, or safety.",
                "Intervals are retrospective cross-validated residual diagnostics, not guaranteed prospective coverage.",
                "Predictions outside the applicability domain must not be treated as reliable.",
            ],
        },
        "source_joblib": {
            "file": source.name,
            "sha256": _sha256(source),
            "conversion_only_no_retraining": True,
        },
    }
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    manifest = {
        "config": config_path.name,
        "config_sha256": _sha256(config_path),
        "model": model_path.name,
        "model_sha256": model_hash,
        "fingerprints": fingerprints_path.name,
        "fingerprints_sha256": fingerprint_hash,
        "model_version": version,
    }
    (output_dir / "inference_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Export a PhytoTruth model to safe native inference artifacts")
    parser.add_argument("--source", type=Path, required=True, help="Trained best_model.joblib")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export_inference_artifacts(args.source, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
