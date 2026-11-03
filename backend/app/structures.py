from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass(frozen=True)
class StandardizedStructure:
    original_smiles: str
    canonical_smiles: str
    standard_inchi: str
    standard_inchi_key: str
    removed_fragments: bool


def standardize_smiles(smiles: str) -> StandardizedStructure:
    """Create a parent, uncharged, stereochemistry-preserving RDKit structure."""
    value = (smiles or "").strip()
    if not value:
        raise ValueError("SMILES is empty")
    molecule = Chem.MolFromSmiles(value, sanitize=True)
    if molecule is None:
        raise ValueError("SMILES could not be parsed by RDKit")
    cleaned = rdMolStandardize.Cleanup(molecule)
    parent = rdMolStandardize.FragmentParent(cleaned)
    parent = rdMolStandardize.Uncharger().uncharge(parent)
    Chem.SanitizeMol(parent)
    canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    if not canonical:
        raise ValueError("RDKit produced an empty parent structure")
    standard_inchi = Chem.MolToInchi(parent)
    inchi_key = Chem.InchiToInchiKey(standard_inchi)
    if not standard_inchi or not inchi_key:
        raise ValueError("RDKit could not generate a standard InChI/InChIKey")
    return StandardizedStructure(
        original_smiles=value,
        canonical_smiles=canonical,
        standard_inchi=standard_inchi,
        standard_inchi_key=inchi_key,
        removed_fragments=len(Chem.GetMolFrags(molecule)) > len(Chem.GetMolFrags(parent)),
    )

