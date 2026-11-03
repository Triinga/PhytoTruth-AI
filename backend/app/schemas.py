from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


class EvidenceRecord(BaseModel):
    molecule_id: str
    activity_id: int | str
    assay_id: str
    document_id: str | None = None
    pchembl_value: float = Field(ge=0, le=15)
    standard_value_nm: float = Field(gt=0)
    relation: Literal["=", "<", ">", "<=", ">="] = "="
    organism: str | None = None
    assay_type: str | None = None
    potential_duplicate: bool = False
    validity_comment: str | None = None


class EvidenceRequest(BaseModel):
    records: list[EvidenceRecord] = Field(min_length=1)
    active_threshold_pchembl: float = 6.0


class EvidenceScore(BaseModel):
    molecule_id: str
    measurements: int
    independent_assays: int
    documents: int
    median_pchembl: float
    max_pchembl: float
    mad: float
    pchembl_range: float
    contradiction_rate: float
    duplicate_rate: float
    organism_count: int
    reliability: float
    evidence_score: float
    ranking_score: float
    ranking_basis: Literal["median_pchembl_development_validated_guardrail"]
    reliability_interpretation: Literal["descriptive_warning_not_probability"]
    confidence_interval: tuple[float, float]
    flags: list[str]


class PredictionRequest(BaseModel):
    model_config = ConfigDict(
        json_schema_extra={"example": {"smiles": "CCO"}},
    )

    smiles: str = Field(
        min_length=1,
        max_length=10_000,
        description="A molecular structure encoded as a SMILES string.",
    )

    @field_validator("smiles")
    @classmethod
    def normalize_input(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("SMILES must not be blank")
        return value


class PredictionInterval(BaseModel):
    lower: float
    upper: float
    confidence_level: float = Field(gt=0, lt=1)
    method: str


class ApplicabilityDomainResult(BaseModel):
    inside: bool
    status: Literal["IN_DOMAIN_RESEARCH_ONLY", "OUTSIDE_DOMAIN_DO_NOT_RELY"]
    nearest_training_tanimoto: float = Field(ge=0, le=1)
    threshold: float = Field(ge=0, le=1)
    warning: str | None = None


class MolecularDescriptors(BaseModel):
    molecular_weight: float
    logp: float
    tpsa: float
    hbd: float
    hba: float
    rotatable_bonds: float
    fraction_csp3: float
    ring_count: float
    heavy_atoms: float
    formal_charge: float


class PredictionModelMetadata(BaseModel):
    name: str
    version: str
    format: str
    target: dict[str, str | None]
    training_structures: int = Field(gt=0)
    training_scaffolds: int = Field(gt=0)


class PredictionResponse(BaseModel):
    canonical_smiles: str
    predicted_pchembl: float
    prediction_interval: PredictionInterval
    applicability_domain: ApplicabilityDomainResult
    descriptors: MolecularDescriptors
    model: PredictionModelMetadata
    kind: Literal["predicted"]
    limitations: list[str]

