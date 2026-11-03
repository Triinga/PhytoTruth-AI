from copy import deepcopy

from app.reliability_investigation import (
    component_predictive_value,
    compound_error_analysis,
    deterministic_partition,
    select_ranking_on_development,
)


def row(molecule_id: str, *, natural: bool = False, median: float = 6.0, full: float = 6.0, truth: float = 6.0):
    return {
        "molecule_id": molecule_id,
        "is_verified_natural_product": natural,
        "train_assays": 2,
        "train_documents": 2,
        "heldout_median_pactivity": truth,
        "median_potency": median,
        "phytotruth_full": full,
        "support_bonus": 0.1,
        "consistency_penalty": 0.2,
        "reliability_penalty": 0.05,
        "uncertainty_penalty": 0.1,
    }


def test_partition_is_deterministic_disjoint_and_stratified():
    rows = [row(f"M{index}", natural=index < 6) for index in range(20)]
    first_dev, first_eval = deterministic_partition(rows, development_fraction=0.6, seed=17)
    second_dev, second_eval = deterministic_partition(list(reversed(rows)), development_fraction=0.6, seed=17)
    assert [item["molecule_id"] for item in first_dev] == [item["molecule_id"] for item in second_dev]
    assert [item["molecule_id"] for item in first_eval] == [item["molecule_id"] for item in second_eval]
    assert {item["molecule_id"] for item in first_dev}.isdisjoint(
        {item["molecule_id"] for item in first_eval}
    )
    assert sum(item["is_verified_natural_product"] for item in first_dev) == 4
    assert sum(item["is_verified_natural_product"] for item in first_eval) == 2


def test_selection_gate_requires_both_metrics_and_ignores_evaluation_data():
    inference = {
        "ndcg_at_k": {"difference_ci95": [0.01, 0.2], "paired_randomization_p": 0.01},
        "kendall_tau": {"difference_ci95": [-0.01, 0.2], "paired_randomization_p": 0.01},
    }
    assert select_ranking_on_development(inference)["selected_method"] == "median_potency"
    inference["kendall_tau"]["difference_ci95"][0] = 0.001
    assert select_ranking_on_development(inference)["selected_method"] == "phytotruth_full"


def test_compound_error_analysis_identifies_harmful_reversal_and_component():
    first = row("A", median=7.0, full=5.7, truth=7.0)
    first.update(consistency_penalty=1.3, support_bonus=0.0, reliability_penalty=0.0, uncertainty_penalty=0.0)
    second = row("B", median=6.0, full=6.0, truth=6.0)
    second.update(consistency_penalty=0.0, support_bonus=0.0, reliability_penalty=0.0, uncertainty_penalty=0.0)
    compounds, pairs, summary = compound_error_analysis([first, second])
    assert summary["harmful_reversals"] == 1
    assert summary["beneficial_reversals"] == 0
    assert pairs[0]["dominant_harmful_component"] == "experimental_consistency"
    assert {item["harmful_pair_count"] for item in compounds} == {1}


def test_component_analysis_is_reproducible_and_heldout_values_are_read_only():
    rows = [
        row(f"M{index}", median=5.0 + index / 10, full=5.0 + index / 10 - index / 100, truth=5.2 + index / 8)
        for index in range(8)
    ]
    original_scores = [(item["median_potency"], item["phytotruth_full"]) for item in rows]
    first = component_predictive_value(rows, iterations=50, seed=9)
    second = component_predictive_value(deepcopy(rows), iterations=50, seed=9)
    assert first == second
    assert original_scores == [(item["median_potency"], item["phytotruth_full"]) for item in rows]
