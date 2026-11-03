from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class RankingMetrics:
    ndcg_at_k: float
    enrichment_factor_at_k: float
    precision_at_k: float
    kendall_tau: float


METRIC_NAMES = tuple(RankingMetrics.__dataclass_fields__)


def _kendall_tau_b(left: np.ndarray, right: np.ndarray) -> float:
    """Kendall's tau-b with explicit tie handling and no optional dependency."""
    if len(left) < 2:
        return 0.0
    first, second = np.triu_indices(len(left), k=1)
    left_sign = np.sign(left[first] - left[second])
    right_sign = np.sign(right[first] - right[second])
    products = left_sign * right_sign
    concordant = int(np.count_nonzero(products > 0))
    discordant = int(np.count_nonzero(products < 0))
    left_only_ties = int(np.count_nonzero((left_sign == 0) & (right_sign != 0)))
    right_only_ties = int(np.count_nonzero((right_sign == 0) & (left_sign != 0)))
    denominator = math.sqrt(
        (concordant + discordant + left_only_ties)
        * (concordant + discordant + right_only_ties)
    )
    return (concordant - discordant) / denominator if denominator else 0.0


def _arrays(scores: Sequence[float], heldout: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    score = np.asarray(scores, dtype=float)
    truth = np.asarray(heldout, dtype=float)
    if score.ndim != 1 or truth.ndim != 1 or len(score) != len(truth):
        raise ValueError("scores and heldout must be one-dimensional arrays of equal length")
    if not len(score):
        raise ValueError("at least one molecule is required")
    if not np.all(np.isfinite(score)) or not np.all(np.isfinite(truth)):
        raise ValueError("scores and heldout values must be finite")
    return score, truth


def evaluate_ranking(
    scores: Sequence[float],
    heldout: Sequence[float],
    k: int = 10,
    *,
    active_threshold: float = 6.0,
) -> RankingMetrics:
    """Evaluate a frozen ranking against evidence that was never used to rank it."""
    score, truth = _arrays(scores, heldout)
    if k < 1:
        raise ValueError("k must be positive")
    k = min(k, len(score))
    # Stable secondary index makes ties reproducible without inspecting held-out truth.
    order = np.lexsort((np.arange(len(score)), -score))
    ideal = np.lexsort((np.arange(len(truth)), -truth))
    # pActivity is logarithmic, so linear gain above the prespecified activity
    # threshold has a direct interpretation in log-potency units.
    gains = np.maximum(truth - active_threshold, 0.0)
    discounts = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float(np.sum(gains[order[:k]] * discounts))
    idcg = float(np.sum(gains[ideal[:k]] * discounts))
    active = truth >= active_threshold
    hit_rate = float(np.mean(active[order[:k]]))
    base_rate = float(np.mean(active))
    tau = _kendall_tau_b(score, truth)
    return RankingMetrics(
        ndcg_at_k=dcg / idcg if idcg else 0.0,
        enrichment_factor_at_k=hit_rate / base_rate if base_rate else 0.0,
        precision_at_k=hit_rate,
        kendall_tau=float(tau),
    )


def compare_rankings(
    potency_scores: Sequence[float],
    evidence_scores: Sequence[float],
    heldout: Sequence[float],
    k: int = 10,
    *,
    active_threshold: float = 6.0,
) -> dict[str, dict[str, float]]:
    return {
        "potency_only": asdict(evaluate_ranking(potency_scores, heldout, k, active_threshold=active_threshold)),
        "evidence_aware": asdict(evaluate_ranking(evidence_scores, heldout, k, active_threshold=active_threshold)),
    }


def _metric(metrics: RankingMetrics, name: str) -> float:
    if name not in METRIC_NAMES:
        raise ValueError(f"Unknown ranking metric: {name}")
    return float(getattr(metrics, name))


def paired_ranking_inference(
    baseline_scores: Sequence[float],
    candidate_scores: Sequence[float],
    heldout: Sequence[float],
    *,
    k: int = 10,
    active_threshold: float = 6.0,
    bootstrap_iterations: int = 2000,
    permutation_iterations: int = 2000,
    seed: int = 2026,
) -> dict[str, dict[str, float | list[float]]]:
    """Molecule-level paired bootstrap CIs and score-swap randomization tests."""
    baseline, truth = _arrays(baseline_scores, heldout)
    candidate, candidate_truth = _arrays(candidate_scores, heldout)
    if not np.array_equal(truth, candidate_truth):
        raise ValueError("heldout arrays differ")
    if bootstrap_iterations < 1 or permutation_iterations < 1:
        raise ValueError("bootstrap and permutation iterations must be positive")

    rng = np.random.default_rng(seed)
    n = len(truth)
    baseline_point = evaluate_ranking(baseline, truth, k, active_threshold=active_threshold)
    candidate_point = evaluate_ranking(candidate, truth, k, active_threshold=active_threshold)
    bootstrap = {name: {"baseline": [], "candidate": [], "difference": []} for name in METRIC_NAMES}
    for _ in range(bootstrap_iterations):
        indices = rng.integers(0, n, size=n)
        left = evaluate_ranking(baseline[indices], truth[indices], k, active_threshold=active_threshold)
        right = evaluate_ranking(candidate[indices], truth[indices], k, active_threshold=active_threshold)
        for name in METRIC_NAMES:
            left_value = _metric(left, name)
            right_value = _metric(right, name)
            bootstrap[name]["baseline"].append(left_value)
            bootstrap[name]["candidate"].append(right_value)
            bootstrap[name]["difference"].append(right_value - left_value)

    observed = {name: _metric(candidate_point, name) - _metric(baseline_point, name) for name in METRIC_NAMES}
    exceedances = {name: 0 for name in METRIC_NAMES}
    for _ in range(permutation_iterations):
        swap = rng.random(n) < 0.5
        permuted_left = np.where(swap, candidate, baseline)
        permuted_right = np.where(swap, baseline, candidate)
        left = evaluate_ranking(permuted_left, truth, k, active_threshold=active_threshold)
        right = evaluate_ranking(permuted_right, truth, k, active_threshold=active_threshold)
        for name in METRIC_NAMES:
            delta = _metric(right, name) - _metric(left, name)
            exceedances[name] += int(abs(delta) >= abs(observed[name]) - 1e-12)

    result: dict[str, dict[str, float | list[float]]] = {}
    for name in METRIC_NAMES:
        samples = bootstrap[name]
        result[name] = {
            "baseline": _metric(baseline_point, name),
            "candidate": _metric(candidate_point, name),
            "difference": observed[name],
            "baseline_ci95": [float(x) for x in np.quantile(samples["baseline"], [0.025, 0.975])],
            "candidate_ci95": [float(x) for x in np.quantile(samples["candidate"], [0.025, 0.975])],
            "difference_ci95": [float(x) for x in np.quantile(samples["difference"], [0.025, 0.975])],
            "paired_randomization_p": (exceedances[name] + 1) / (permutation_iterations + 1),
        }
    return result
