import numpy as np
import pytest

from app.features import DESCRIPTOR_NAMES, feature_config, molecular_feature_vector


def test_feature_vector_contains_fingerprint_and_ordered_descriptors():
    vector = molecular_feature_vector("CCO", radius=2, n_bits=128)
    config = feature_config(radius=2, n_bits=128)
    assert vector.shape == (128 + len(DESCRIPTOR_NAMES),)
    assert vector.dtype == np.float32
    assert np.all(np.isfinite(vector))
    assert config["descriptor_names"] == list(DESCRIPTOR_NAMES)
    assert config["feature_count"] == len(vector)


def test_feature_generation_is_deterministic_and_validates_smiles():
    assert np.array_equal(molecular_feature_vector("CCO"), molecular_feature_vector("OCC"))
    with pytest.raises(ValueError, match="Invalid SMILES"):
        molecular_feature_vector("not-a-smiles")

