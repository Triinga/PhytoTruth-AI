import math

from app.robust_multitarget_evaluation import (
    _exact_sign_flip_p,
    _holm_adjust,
    deterministic_target_partition,
    target_cluster_inference,
)


def _row(molecule, median, candidate, truth):
    return {
        "molecule_id": molecule,
        "median_potency": median,
        "optimized_ranking": candidate,
        "heldout_median_pactivity": truth,
    }


def test_target_partition_is_deterministic_and_disjoint():
    targets = [f"T{index}" for index in range(20)]
    first = deterministic_target_partition(targets, development_fraction=0.6, seed=2026)
    second = deterministic_target_partition(targets, development_fraction=0.6, seed=2026)
    assert first == second
    assert set(first[0]).isdisjoint(first[1])
    assert set(first[0]) | set(first[1]) == set(targets)


def test_exact_sign_flip_has_target_level_resolution():
    differences = [0.1] * 11
    assert _exact_sign_flip_p(differences) == 2 / (2 ** 11)


def test_holm_adjustment_is_monotone_and_bounded():
    adjusted = _holm_adjust({"a": 0.01, "b": 0.03, "c": 0.2})
    assert adjusted == {"a": 0.03, "b": 0.06, "c": 0.2}


def test_target_cluster_inference_does_not_count_molecules_as_targets():
    rows = {
        "T1": [_row("A", 7, 7, 7), _row("B", 6, 6, 6), _row("C", 5, 5, 5)],
        "T2": [_row("D", 5, 7, 7), _row("E", 7, 5, 5), _row("F", 6, 6, 6)],
    }
    result = target_cluster_inference(
        rows,
        "optimized_ranking",
        k=2,
        active_threshold=6,
        bootstrap_iterations=100,
        seed=7,
    )
    assert result["kendall_tau"]["difference"] > 0
    assert result["kendall_tau"]["target_wins"] == 1
    assert result["kendall_tau"]["target_ties"] == 1
    assert math.isclose(result["kendall_tau"]["target_sign_flip_p"], 1.0)
