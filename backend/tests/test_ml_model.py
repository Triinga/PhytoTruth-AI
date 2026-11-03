from pathlib import Path

import joblib
import numpy as np
import pytest

from app.features import feature_config, molecular_feature_vector
from app.ml import (
    load_model,
    make_regressor,
    nearest_neighbor_threshold,
    pack_fingerprints,
    predict_activity,
)
from app.ml_experiment import assign_scaffold_folds


SMILES = ["c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "C1CCCCC1", "CCO", "CCN"]
TARGETS = np.asarray([5.0, 5.4, 5.8, 6.1, 4.7, 4.9])


def test_scaffold_assignment_is_reproducible_and_has_no_group_leakage():
    scaffolds = ["benzene", "benzene", "benzene", "cyclohexane", "ethanol", "ethylamine"]
    first = assign_scaffold_folds(scaffolds, folds=3)
    second = assign_scaffold_folds(scaffolds, folds=3)
    assert np.array_equal(first, second)
    for scaffold in set(scaffolds):
        assert len(set(first[index] for index, value in enumerate(scaffolds) if value == scaffold)) == 1


def test_random_forest_training_is_reproducible():
    X = np.stack([molecular_feature_vector(smiles, n_bits=128) for smiles in SMILES])
    first = make_regressor("random_forest", seed=19)
    second = make_regressor("random_forest", seed=19)
    first.fit(X, TARGETS)
    second.fit(X, TARGETS)
    assert np.array_equal(first.predict(X), second.predict(X))


def test_serialized_model_loads_and_predicts_with_uncertainty_and_domain(tmp_path: Path):
    config = feature_config(n_bits=128)
    X = np.stack([molecular_feature_vector(smiles, n_bits=128) for smiles in SMILES])
    model = make_regressor("random_forest", seed=23)
    model.fit(X, TARGETS)
    packed = pack_fingerprints(X, n_bits=128)
    bundle = {
        "model": model,
        "feature_config": config,
        "uncertainty": {"conformal_half_width_90": 0.8},
        "applicability_domain": {
            "similarity_threshold": nearest_neighbor_threshold(packed),
            "reference_packed_fingerprints": packed,
        },
        "training": {"model_name": "random_forest"},
    }
    path = tmp_path / "model.joblib"
    joblib.dump(bundle, path)
    loaded = load_model(path)
    result = predict_activity(loaded, "CCO")
    assert result["kind"] == "predicted"
    assert result["prediction_interval_90"][0] < result["predicted_pactivity"] < result["prediction_interval_90"][1]
    assert 0 <= result["nearest_training_tanimoto"] <= 1
    assert isinstance(result["inside_applicability_domain"], bool)


def test_invalid_model_bundle_is_rejected(tmp_path: Path):
    path = tmp_path / "bad.joblib"
    joblib.dump({"model": None}, path)
    with pytest.raises(ValueError, match="Invalid PhytoTruth model bundle"):
        load_model(path)
