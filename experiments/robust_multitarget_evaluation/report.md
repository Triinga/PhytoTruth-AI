# Robust multi-target PhytoTruth evaluation

Generated 2026-10-26T14:28:55.200516+00:00 from real ChEMBL IC50 measurements already present in the repository.

## Leakage contract

Whole biological targets were deterministically assigned to development or evaluation before any outcome-based model fitting. For every molecule, the latest publication year was held out and all scoring features came from strictly earlier, document-disjoint evidence. Ridge strength and feature family were selected only by leave-one-development-target-out validation. The evaluation targets were used once for the frozen comparison.

This is a retrospective database benchmark, not a prospective or independent-laboratory replication. Several source targets have appeared in earlier repository analyses, so the result is not described as a new confirmatory study.

## Coverage

- Discovered targets: 30
- Eligible targets: 20
- Development targets: 12
- Evaluation targets: 8
- Evaluation molecules with future-publication evidence: 1235

## Development-only selection

The best nested-development candidate was `reliability:ridge:10`. It did not pass the frozen two-metric gate, so the deployed ranking is `median_potency`.

## Locked evaluation

Optimized candidate minus median potency macro nDCG@10: 0.0007 (95% target-bootstrap CI -0.0053 to 0.0084; exact target sign-flip p=0.890625; Holm-adjusted p=1.000000; Hedges paired g=0.060; wins/ties/losses=2/0/6).

Kendall tau difference: -0.0088 (95% CI -0.0160 to -0.0014; exact p=0.062500; Holm-adjusted p=0.250000).

Scientific status: **superiority_not_demonstrated**. The frozen target-held-out evaluation does not support a superiority claim; the median-potency guardrail remains deployed.

## Statistical power

With 8 target clusters, the minimum attainable two-sided exact sign-flip p-value is 0.007812. The observed paired standardized effect was 0.060. Approximate target count for 80% power at alpha 0.05: 2157.

## Ablations

`ablation_inference.json` and `ablation_summary.csv` compare the frozen full candidate with removal of support, consistency, reliability/uncertainty, and recent-drift features. These are secondary analyses; Holm-adjusted target-level p-values are reported and negative contributions are retained.

## Interpretation limits

- ChEMBL publications are proxies for independent experiments and may share laboratories or copied measurements.
- Target-level inference avoids treating molecules as independent biological targets, but available targets are still a convenience sample.
- nDCG@10 measures future-publication ranking agreement, not efficacy, safety, or clinical utility.
- Natural-product membership is never a training feature; subgroup counts are descriptive because target-held-out NP power remains limited.
