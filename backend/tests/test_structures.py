import pytest

from app.structures import standardize_smiles


def test_salt_parent_and_neutral_form_deduplicate():
    sodium_acetate = standardize_smiles("CC(=O)[O-].[Na+]")
    acetic_acid = standardize_smiles("CC(=O)O")
    assert sodium_acetate.standard_inchi_key == acetic_acid.standard_inchi_key
    assert sodium_acetate.removed_fragments is True


def test_stereochemistry_is_preserved():
    left = standardize_smiles("N[C@@H](C)C(=O)O")
    right = standardize_smiles("N[C@H](C)C(=O)O")
    assert left.standard_inchi_key != right.standard_inchi_key
    assert "@" in left.canonical_smiles


@pytest.mark.parametrize("value", ["", "not-a-smiles", "C1CC"])
def test_invalid_smiles_are_rejected(value):
    with pytest.raises(ValueError):
        standardize_smiles(value)

