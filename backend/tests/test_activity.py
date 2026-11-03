import pytest

from app.activity import normalize_chembl_activity


def activity(**overrides):
    record = {
        "activity_id": 101,
        "molecule_chembl_id": "CHEMBL1",
        "canonical_smiles": "CCO",
        "assay_chembl_id": "CHEMBL_A1",
        "document_chembl_id": "CHEMBL_D1",
        "target_chembl_id": "CHEMBL335",
        "target_organism": "Homo sapiens",
        "standard_type": "IC50",
        "standard_relation": "=",
        "standard_value": "100",
        "standard_units": "nM",
        "pchembl_value": "7.0",
        "assay_type": "B",
        "bao_format": "BAO_0000357",
    }
    record.update(overrides)
    return record


def test_units_are_normalized_to_molar_and_pactivity():
    result = normalize_chembl_activity(activity()).record
    assert result is not None
    assert result["standard_value_molar"] == pytest.approx(1e-7)
    assert result["computed_pactivity"] == pytest.approx(7.0)
    assert result["pchembl_consistent"] is True


def test_incompatible_endpoint_or_format_gets_distinct_context():
    ic50 = normalize_chembl_activity(activity()).record
    ki = normalize_chembl_activity(activity(standard_type="Ki")).record
    cellular = normalize_chembl_activity(activity(bao_format="BAO_0000219")).record
    assert ic50 and ki and cellular
    assert len({ic50["assay_context_key"], ki["assay_context_key"], cellular["assay_context_key"]}) == 3


def test_censored_value_is_preserved_not_treated_as_exact():
    result = normalize_chembl_activity(activity(standard_relation=">", pchembl_value="7.0")).record
    assert result is not None
    assert result["is_censored"] is True
    assert result["pchembl_consistent"] is None


def test_non_concentration_units_are_rejected():
    result = normalize_chembl_activity(activity(standard_units="%"))
    assert result.record is None
    assert "unsupported" in result.rejection_reason

