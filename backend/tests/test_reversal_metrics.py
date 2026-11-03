import pytest

from app.reversal import evaluate_ranking, paired_ranking_inference


def test_perfect_ranking_has_unit_ndcg_and_precision():
    metrics = evaluate_ranking([3, 2, 1], [8, 7, 5], k=2, active_threshold=6)
    assert metrics.ndcg_at_k == pytest.approx(1.0)
    assert metrics.precision_at_k == pytest.approx(1.0)
    assert metrics.kendall_tau == pytest.approx(1.0)


def test_ranking_validation_rejects_bad_inputs():
    with pytest.raises(ValueError, match="equal length"):
        evaluate_ranking([1], [1, 2])
    with pytest.raises(ValueError, match="positive"):
        evaluate_ranking([1], [1], k=0)


def test_paired_inference_is_reproducible():
    kwargs = dict(
        baseline_scores=[3, 2, 1, 0],
        candidate_scores=[4, 1, 3, 0],
        heldout=[8, 5, 7, 4],
        k=2,
        bootstrap_iterations=30,
        permutation_iterations=30,
        seed=17,
    )
    assert paired_ranking_inference(**kwargs) == paired_ranking_inference(**kwargs)

