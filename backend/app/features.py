from __future__ import annotations

import numpy as np

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, Lipinski
    from rdkit.Chem.Scaffolds import MurckoScaffold
except ImportError as exc:  # pragma: no cover - installation check
    raise RuntimeError("RDKit is required. Install backend/requirements.txt.") from exc


DESCRIPTOR_NAMES = (
    "molecular_weight",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "fraction_csp3",
    "ring_count",
    "heavy_atoms",
    "formal_charge",
)


def canonicalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    generator = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = generator.GetFingerprint(mol)
    array = np.zeros((n_bits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, array)
    return array


def molecular_descriptors(smiles: str) -> dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    return {
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "logp": round(Descriptors.MolLogP(mol), 3),
        "tpsa": round(Descriptors.TPSA(mol), 3),
        "hbd": float(Lipinski.NumHDonors(mol)),
        "hba": float(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
        "fraction_csp3": float(Lipinski.FractionCSP3(mol)),
        "ring_count": float(Lipinski.RingCount(mol)),
        "heavy_atoms": float(Lipinski.HeavyAtomCount(mol)),
        "formal_charge": float(Chem.GetFormalCharge(mol)),
    }


def molecular_feature_vector(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Morgan bits followed by a fixed, ordered RDKit descriptor vector."""
    fingerprint = morgan_fingerprint(smiles, radius=radius, n_bits=n_bits).astype(np.float32)
    descriptors = molecular_descriptors(smiles)
    descriptor_vector = np.asarray([descriptors[name] for name in DESCRIPTOR_NAMES], dtype=np.float32)
    if not np.all(np.isfinite(descriptor_vector)):
        raise ValueError("RDKit produced a non-finite molecular descriptor")
    return np.concatenate((fingerprint, descriptor_vector))


def feature_config(radius: int = 2, n_bits: int = 2048) -> dict[str, object]:
    return {
        "fingerprint": "Morgan",
        "radius": radius,
        "n_bits": n_bits,
        "use_chirality": False,
        "descriptor_names": list(DESCRIPTOR_NAMES),
        "feature_count": n_bits + len(DESCRIPTOR_NAMES),
    }


def bemis_murcko_scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    value = Chem.MolToSmiles(scaffold)
    return value or canonicalize_smiles(smiles)

