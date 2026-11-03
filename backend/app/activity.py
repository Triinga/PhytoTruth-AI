from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any


UNIT_TO_MOLAR = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "nM": 1e-9, "pM": 1e-12}
EXACT_RELATIONS = {"=", None, ""}


@dataclass(frozen=True)
class NormalizationResult:
    record: dict[str, Any] | None
    rejection_reason: str | None = None


def assay_context(record: dict[str, Any]) -> dict[str, str | None]:
    """Fields that must match before quantitative measurements are aggregated."""
    return {
        "target_chembl_id": record.get("target_chembl_id"),
        "target_organism": record.get("target_organism"),
        "standard_type": record.get("standard_type"),
        "assay_type": record.get("assay_type"),
        "bao_format": record.get("bao_format"),
    }


def context_key(context: dict[str, Any]) -> str:
    encoded = json.dumps(context, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(encoded.encode()).hexdigest()[:20]


def normalize_chembl_activity(record: dict[str, Any]) -> NormalizationResult:
    """Normalize molar concentration endpoints without merging incompatible contexts."""
    activity_id = record.get("activity_id")
    if activity_id in (None, ""):
        return NormalizationResult(None, "missing activity_id")
    unit = record.get("standard_units")
    if unit not in UNIT_TO_MOLAR:
        return NormalizationResult(None, f"unsupported standard_units: {unit!r}")
    try:
        value = float(record.get("standard_value"))
    except (TypeError, ValueError):
        return NormalizationResult(None, "standard_value is not numeric")
    if not math.isfinite(value) or value <= 0:
        return NormalizationResult(None, "standard_value must be finite and positive")
    standard_type = str(record.get("standard_type") or "").strip()
    target = str(record.get("target_chembl_id") or "").strip()
    assay_id = str(record.get("assay_chembl_id") or "").strip()
    if not standard_type or not target or not assay_id:
        return NormalizationResult(None, "missing endpoint, target, or assay identifier")

    molar = value * UNIT_TO_MOLAR[unit]
    computed_pactivity = -math.log10(molar)
    reported_pchembl = record.get("pchembl_value")
    reported = None
    if reported_pchembl not in (None, ""):
        try:
            reported = float(reported_pchembl)
        except (TypeError, ValueError):
            return NormalizationResult(None, "pchembl_value is not numeric")
    relation = record.get("standard_relation") or record.get("relation") or "="
    context = assay_context(record)
    normalized = {
        "source": "ChEMBL",
        "source_record_id": str(activity_id),
        "activity_id": activity_id,
        "molecule_chembl_id": record.get("molecule_chembl_id"),
        "canonical_smiles": record.get("canonical_smiles"),
        "assay_chembl_id": assay_id,
        "document_chembl_id": record.get("document_chembl_id"),
        "target_chembl_id": target,
        "target_organism": record.get("target_organism"),
        "standard_type": standard_type,
        "standard_relation": relation,
        "original_standard_value": value,
        "original_standard_units": unit,
        "standard_value_molar": molar,
        "computed_pactivity": computed_pactivity,
        "reported_pchembl_value": reported,
        "pchembl_consistent": None if reported is None or relation not in EXACT_RELATIONS else abs(reported - computed_pactivity) <= 0.15,
        "is_censored": relation not in EXACT_RELATIONS,
        "assay_context": context,
        "assay_context_key": context_key(context),
        "data_validity_comment": record.get("data_validity_comment"),
        "potential_duplicate": bool(record.get("potential_duplicate")),
        "raw_source": record,
    }
    return NormalizationResult(normalized)

