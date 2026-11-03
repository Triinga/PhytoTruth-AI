from __future__ import annotations

import hashlib
import math
from collections import defaultdict

import numpy as np

from .schemas import EvidenceRecord, EvidenceScore


def _stable_seed(value: str) -> int:
    return int(hashlib.sha256(value.encode()).hexdigest()[:8], 16)


def _bootstrap_interval(values: np.ndarray, molecule_id: str) -> tuple[float, float]:
    if len(values) == 1:
        return (round(float(values[0]), 3), round(float(values[0]), 3))
    rng = np.random.default_rng(_stable_seed(molecule_id))
    samples = rng.choice(values, size=(1000, len(values)), replace=True)
    medians = np.median(samples, axis=1)
    low, high = np.quantile(medians, [0.025, 0.975])
    return round(float(low), 3), round(float(high), 3)


def score_evidence(records: list[EvidenceRecord], threshold: float = 6.0) -> list[EvidenceScore]:
    """Aggregate experimental evidence without silently treating replicates as truth.

    Reliability rewards independent assays/documents and penalises robust dispersion,
    threshold-crossing contradictions, duplicates, invalidity comments and organism mixing.
    The legacy evidence score is retained for reproducibility. Ranking uses the robust
    median guardrail because the reliability investigation found that dispersion and
    support predict error magnitude, but not the direction needed for a potency
    adjustment. Reliability is descriptive and is not a calibrated probability.
    """
    grouped: dict[str, list[EvidenceRecord]] = defaultdict(list)
    for record in records:
        grouped[record.molecule_id].append(record)

    results: list[EvidenceScore] = []
    for molecule_id, group in grouped.items():
        values = np.asarray([r.pchembl_value for r in group], dtype=float)
        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)))
        spread = float(np.ptp(values)) if len(values) > 1 else 0.0
        actives = values >= threshold
        contradictory_pairs = sum(
            bool(actives[i] != actives[j])
            for i in range(len(values))
            for j in range(i + 1, len(values))
        )
        pairs = len(values) * (len(values) - 1) / 2
        contradiction = contradictory_pairs / pairs if pairs else 0.0
        assays = {r.assay_id for r in group}
        documents = {r.document_id for r in group if r.document_id}
        organisms = {r.organism for r in group if r.organism}
        duplicate_rate = sum(r.potential_duplicate for r in group) / len(group)
        invalid_rate = sum(bool(r.validity_comment) for r in group) / len(group)

        support = 1 - math.exp(-len(assays) / 2)
        consistency = math.exp(-0.9 * mad) * (1 - 0.55 * contradiction)
        provenance = (0.65 + 0.35 * (1 - math.exp(-len(documents)))) if documents else 0.6
        organism_penalty = 1 / (1 + 0.12 * max(0, len(organisms) - 1))
        reliability = support * consistency * provenance * organism_penalty
        reliability *= 1 - 0.45 * duplicate_rate
        reliability *= 1 - 0.55 * invalid_rate
        reliability = float(np.clip(reliability, 0, 1))

        # Conservative lower-bound utility: central potency minus inconsistency and
        # uncertainty penalties, plus a modest independent-assay support reward.
        evidence_score = (
            median
            - 0.65 * mad
            - 0.8 * contradiction
            - 0.12 * spread
            + 0.18 * math.log1p(len(assays))
        )
        flags: list[str] = []
        if spread >= 1.5:
            flags.append("high assay spread")
        if contradiction > 0:
            flags.append("activity threshold reversal")
        if len(organisms) > 1:
            flags.append("cross-organism evidence")
        if duplicate_rate:
            flags.append("potential duplicate")
        if invalid_rate:
            flags.append("source validity warning")
        if len(assays) < 2:
            flags.append("single-assay evidence")

        results.append(EvidenceScore(
            molecule_id=molecule_id,
            measurements=len(group),
            independent_assays=len(assays),
            documents=len(documents),
            median_pchembl=round(median, 3),
            max_pchembl=round(float(np.max(values)), 3),
            mad=round(mad, 3),
            pchembl_range=round(spread, 3),
            contradiction_rate=round(contradiction, 3),
            duplicate_rate=round(duplicate_rate, 3),
            organism_count=len(organisms),
            reliability=round(reliability, 3),
            evidence_score=round(evidence_score, 3),
            ranking_score=round(median, 3),
            ranking_basis="median_pchembl_development_validated_guardrail",
            reliability_interpretation="descriptive_warning_not_probability",
            confidence_interval=_bootstrap_interval(values, molecule_id),
            flags=flags,
        ))
    return sorted(results, key=lambda item: (item.ranking_score, item.molecule_id), reverse=True)

