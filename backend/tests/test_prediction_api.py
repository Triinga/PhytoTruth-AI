from __future__ import annotations

import json
from pathlib import Path

import pytest
from fastapi import HTTPException
from pydantic import ValidationError

from app.main import app, predict_smiles
from app.prediction_service import clear_model_cache, load_predictor
from app.schemas import PredictionRequest


MODEL_DIRECTORY = Path(__file__).resolve().parents[2] / "experiments" / "ml_activity_chembl335"


@pytest.fixture(autouse=True)
def use_test_model(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("PHYTOTRUTH_MODEL_DIR", str(MODEL_DIRECTORY))
    clear_model_cache()
    yield
    clear_model_cache()


def predict(smiles: str) -> dict:
    return predict_smiles(PredictionRequest(smiles=smiles))


def test_predicts_with_real_native_model_and_reports_uncertainty():
    result = predict(" CCO ")
    assert result["canonical_smiles"] == "CCO"
    assert isinstance(result["predicted_pchembl"], float)
    assert result["prediction_interval"]["lower"] < result["predicted_pchembl"]
    assert result["prediction_interval"]["upper"] > result["predicted_pchembl"]
    assert result["prediction_interval"]["confidence_level"] == pytest.approx(0.9)
    assert result["applicability_domain"]["status"] in {
        "IN_DOMAIN_RESEARCH_ONLY",
        "OUTSIDE_DOMAIN_DO_NOT_RELY",
    }
    assert 0 <= result["applicability_domain"]["nearest_training_tanimoto"] <= 1
    assert result["descriptors"]["molecular_weight"] == pytest.approx(46.069, abs=0.001)
    assert result["model"]["name"] == "xgboost"
    assert result["model"]["format"] == "xgboost-json"
    assert result["model"]["target"]["target_chembl_id"] == "CHEMBL335"
    assert result["kind"] == "predicted"
    assert any("clinical" in limitation.lower() for limitation in result["limitations"])


def test_prediction_is_deterministic_and_model_is_cached():
    first = predict("CCOc1ccccc1")
    second = predict("CCOc1ccccc1")

    assert first == second
    cache = load_predictor.cache_info()
    assert cache.misses == 1
    assert cache.hits >= 1


@pytest.mark.parametrize("smiles", ["not-a-smiles", "C1CC"])
def test_invalid_smiles_returns_meaningful_422(smiles: str):
    with pytest.raises(HTTPException) as error:
        predict(smiles)

    assert error.value.status_code == 422
    assert error.value.detail["code"] == "invalid_smiles"


def test_blank_smiles_is_rejected_by_request_validation():
    with pytest.raises(ValidationError):
        PredictionRequest(smiles="   ")


def test_missing_model_returns_service_unavailable(
    monkeypatch: pytest.MonkeyPatch,
):
    missing_directory = Path(__file__).resolve().parent / "definitely-missing-model-directory"
    assert not missing_directory.exists()
    monkeypatch.setenv("PHYTOTRUTH_MODEL_DIR", str(missing_directory))
    clear_model_cache()

    with pytest.raises(HTTPException) as error:
        predict("CCO")

    assert error.value.status_code == 503
    detail = error.value.detail
    assert detail["code"] == "model_unavailable"
    assert "Missing model configuration" in detail["message"]


def test_native_artifact_metadata_and_checksums_are_loaded():
    predictor = load_predictor(str(MODEL_DIRECTORY.resolve()))
    config = predictor.config

    assert config["artifact_version"] == 2
    assert config["model"]["format"] == "xgboost-json"
    assert config["source_joblib"]["conversion_only_no_retraining"] is True
    assert predictor.model.n_features_in_ == config["feature_config"]["feature_count"]
    with (MODEL_DIRECTORY / "inference_manifest.json").open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    assert manifest["model_version"] == config["model"]["version"]


def test_outside_domain_prediction_is_explicitly_unreliable():
    result = predict("[He]")
    assert result["applicability_domain"]["inside"] is False
    assert result["applicability_domain"]["status"] == "OUTSIDE_DOMAIN_DO_NOT_RELY"
    assert result["applicability_domain"]["warning"]
    assert result["limitations"][0].startswith("OUTSIDE APPLICABILITY DOMAIN")


def test_predict_route_is_published_in_openapi_schema():
    operation = app.openapi()["paths"]["/api/predict"]["post"]

    assert operation["requestBody"]["required"] is True
    assert operation["responses"]["200"]["content"]["application/json"]["schema"]
    assert "research estimate" in operation["description"].lower()
