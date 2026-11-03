from app.evidence import score_evidence
from app.schemas import EvidenceRecord


def record(activity_id: int, assay: str, value: float, duplicate: bool = False):
    return EvidenceRecord(
        molecule_id="M1",
        activity_id=activity_id,
        assay_id=assay,
        document_id=f"D{activity_id}",
        pchembl_value=value,
        standard_value_nm=10 ** (9 - value),
        organism="Homo sapiens",
        potential_duplicate=duplicate,
    )


def test_flags_threshold_reversal_and_spread():
    result = score_evidence([record(1, "A1", 5.5), record(2, "A2", 8.5)])[0]
    assert "activity threshold reversal" in result.flags
    assert "high assay spread" in result.flags
    assert result.reliability < 0.6


def test_consistent_independent_evidence_scores_higher():
    consistent = score_evidence([record(1, "A1", 7.2), record(2, "A2", 7.3), record(3, "A3", 7.1)])[0]
    noisy = score_evidence([record(1, "A1", 5.2), record(2, "A2", 8.8), record(3, "A3", 7.2)])[0]
    assert consistent.reliability > noisy.reliability
    assert consistent.evidence_score > noisy.evidence_score


def test_validated_ranking_separates_central_potency_from_reliability_warning():
    result = score_evidence([record(1, "A1", 5.2), record(2, "A2", 8.8), record(3, "A3", 7.2)])[0]
    assert result.ranking_score == result.median_pchembl == 7.2
    assert result.ranking_basis == "median_pchembl_development_validated_guardrail"
    assert result.reliability_interpretation == "descriptive_warning_not_probability"


def test_duplicate_record_is_penalized():
    base = score_evidence([record(1, "A1", 7.0), record(2, "A2", 7.0)])[0]
    duplicate = score_evidence([record(1, "A1", 7.0), record(2, "A2", 7.0, True)])[0]
    assert duplicate.reliability < base.reliability

