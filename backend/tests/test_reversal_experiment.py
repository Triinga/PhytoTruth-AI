from copy import deepcopy

import pytest

from app.reversal_experiment import (
    build_benchmark_rows,
    global_publication_split,
    latest_publication_split,
    select_compatible_measurements,
)


def measurement(activity_id, molecule, year, value, *, context="CTX", document=None, assay=None):
    return {
        "activity_id": activity_id,
        "source_record_id": str(activity_id),
        "molecule_chembl_id": molecule,
        "computed_pactivity": value,
        "standard_value_molar": 10 ** (-value),
        "standard_relation": "=",
        "is_censored": False,
        "assay_chembl_id": assay or f"A{activity_id}",
        "document_chembl_id": document or f"D{activity_id}",
        "target_organism": "Homo sapiens",
        "potential_duplicate": False,
        "data_validity_comment": None,
        "assay_context_key": context,
        "assay_context": {
            "target_chembl_id": "CHEMBL335",
            "target_organism": "Homo sapiens",
            "standard_type": "IC50",
            "assay_type": "B",
            "bao_format": "BAO_0000357",
        },
        "raw_source": {"document_year": year},
        "document_year": year,
    }


def test_latest_split_is_strict_and_document_disjoint():
    records = [
        measurement(1, "M1", 2010, 6.0),
        measurement(2, "M1", 2012, 7.0),
        measurement(3, "M1", 2012, 7.2),
        measurement(4, "M2", 2011, 5.0),
    ]
    split, exclusions = latest_publication_split(records)
    assert set(split) == {"M1"}
    assert {row["activity_id"] for row in split["M1"]["train"]} == {1}
    assert {row["activity_id"] for row in split["M1"]["heldout"]} == {2, 3}
    assert exclusions["single_publication_year"] == 1


def test_global_split_never_reuses_a_document():
    records = [measurement(1, "M1", 2010, 6.0, document="D1"), measurement(2, "M1", 2020, 7.0, document="D2")]
    split, _ = global_publication_split(records, 2015)
    assert split["M1"]["train"][0]["document_chembl_id"] == "D1"
    assert split["M1"]["heldout"][0]["document_chembl_id"] == "D2"


def test_heldout_values_cannot_change_ranking_scores():
    records = [
        measurement(1, "M1", 2010, 6.1),
        measurement(2, "M1", 2011, 6.4),
        measurement(3, "M1", 2020, 7.5),
    ]
    split, _ = latest_publication_split(records)
    original = build_benchmark_rows(split, {})[0]
    changed_split = deepcopy(split)
    changed_split["M1"]["heldout"][0]["computed_pactivity"] = 3.0
    changed = build_benchmark_rows(changed_split, {})[0]
    for method in (
        "baseline_max_potency",
        "median_potency",
        "mean_potency",
        "latest_training_potency",
        "lower_confidence_bound",
        "ablate_reliability",
        "ablate_uncertainty",
        "ablate_explicit_consistency",
        "ablate_support",
        "phytotruth_full",
    ):
        assert changed[method] == original[method]
    assert changed["heldout_median_pactivity"] != original["heldout_median_pactivity"]


def test_context_selection_uses_coverage_not_activity_values():
    records = [
        measurement(1, "M1", 2010, 2.0, context="BIG"),
        measurement(2, "M2", 2010, 2.0, context="BIG"),
        measurement(3, "M3", 2010, 12.0, context="SMALL"),
    ]
    selected, context, exclusions = select_compatible_measurements(records)
    assert context == "BIG"
    assert {row["molecule_chembl_id"] for row in selected} == {"M1", "M2"}
    assert exclusions["other_assay_context"] == 1


def test_split_rejects_document_leakage():
    records = [
        measurement(1, "M1", 2010, 6.0, document="SAME"),
        measurement(2, "M1", 2020, 7.0, document="SAME"),
    ]
    with pytest.raises(ValueError, match="Document leakage"):
        latest_publication_split(records)
