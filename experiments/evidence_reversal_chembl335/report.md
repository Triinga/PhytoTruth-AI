# PhytoTruth Evidence Reversal Experiment

Generated: 2026-10-20T11:55:14.956329+00:00

## Research question

Does a prespecified evidence-aware heuristic rank compounds by potency that reproduces in later publications better than ranking by the strongest earlier measurement?

## Data and protocol

- Source: normalized real ChEMBL CHEMBL335 IC50 records; NPASS membership is used for the natural-product subgroup.
- Compatible context: `bee3486ec6fbd9a6ba15` ({"assay_type": "B", "bao_format": "BAO_0000357", "standard_type": "IC50", "target_chembl_id": "CHEMBL335", "target_organism": "Homo sapiens"}).
- Only exact, dated measurements were eligible. Censored observations and other assay contexts were excluded rather than pooled.
- Primary split: for each molecule, all measurements from its latest publication year are held out; only strictly earlier documents score the molecule.
- Sensitivity split: global publication cutoff 2015 (training through the cutoff, later evidence held out).
- Activity threshold: pActivity ≥ 6.0; K=10.
- Inference: 2000 molecule-level paired bootstrap samples and 2000 paired score-swap randomizations.
- The reliability value is a deterministic heuristic, not a calibrated probability. No weight or threshold was tuned using held-out evidence.

## Primary results

| Cohort | Molecules | Method | nDCG@K | EF@K | Precision@K | Kendall τ |
|---|---:|---|---:|---:|---:|---:|
| latest_all | 162 | baseline_max_potency | 0.891 | 1.519 | 0.900 | 0.845 |
| latest_all | 162 | phytotruth_full | 0.977 | 1.688 | 1.000 | 0.867 |
| latest_natural | 31 | baseline_max_potency | 0.891 | 2.325 | 0.600 | 0.527 |
| latest_natural | 31 | phytotruth_full | 0.765 | 1.938 | 0.500 | 0.510 |
| global_all | 39 | baseline_max_potency | 0.835 | 1.733 | 0.800 | 0.613 |
| global_all | 39 | phytotruth_full | 0.985 | 2.167 | 1.000 | 0.653 |
| global_natural | 14 | baseline_max_potency | 0.647 | 1.400 | 0.400 | 0.363 |
| global_natural | 14 | phytotruth_full | 0.863 | 1.050 | 0.300 | 0.341 |

For the primary all-compound cohort, the PhytoTruth-minus-baseline nDCG difference was 0.086 (95% molecule-bootstrap CI [0.000, 0.265], paired randomization p=0.3458).

The evidence-aware heuristic's all-compound nDCG point estimate was higher, but its paired 95% interval included zero and the randomization test was not significant. The experiment therefore does not demonstrate an advantage over potency-only ranking.

In the NPASS-verified natural-product cohort, PhytoTruth-minus-baseline nDCG was -0.126 (95% CI [-0.758, 0.418]); the available natural-product data do not support an improvement.

## Ablation study

| Method | nDCG@K | EF@K | Precision@K | Kendall τ |
|---|---:|---:|---:|---:|
| baseline_max_potency | 0.891 | 1.519 | 0.900 | 0.845 |
| median_potency | 0.977 | 1.688 | 1.000 | 0.901 |
| ablate_reliability | 0.977 | 1.688 | 1.000 | 0.868 |
| ablate_uncertainty | 0.977 | 1.688 | 1.000 | 0.871 |
| ablate_explicit_consistency | 0.977 | 1.688 | 1.000 | 0.883 |
| ablate_support | 0.977 | 1.688 | 1.000 | 0.869 |
| phytotruth_full | 0.977 | 1.688 | 1.000 | 0.867 |

Ablations remove one prespecified component while keeping the split fixed. They are diagnostic comparisons, not independently tuned models.

Median potency alone matched the full heuristic's top-K metrics and its Kendall τ exceeded the full heuristic by 0.035. The full-minus-median Kendall difference was -0.035 (95% CI [-0.079, -0.003], p=0.0400). Thus the reliability, uncertainty, consistency, and support terms did not improve this benchmark beyond robust central potency and may have modestly worsened global ordering.

## Leakage audit

- Activity IDs and document IDs are disjoint between scoring and held-out evidence for every molecule.
- Every held-out publication year is strictly later than that molecule's training years.
- Assay context selection used molecule coverage only, never activity outcomes.
- Confidence intervals resample whole molecules; individual measurements are never treated as independent experimental units.

## Limitations

- Only 31 NPASS-verified natural products had repeated evidence across publication years in the selected context. Natural-product subgroup estimates are therefore exploratory and may be very imprecise.
- The latest-year split is retrospective and molecule-specific; it is not a prospective trial. The global-cutoff sensitivity is stricter but smaller.
- ChEMBL assay-context metadata cannot capture every protocol difference. Sharing target, organism, endpoint, assay type, and BAO format reduces but does not eliminate heterogeneity.
- Held-out ChEMBL measurements are literature evidence, not clinical efficacy or safety outcomes.
- The heuristic weights were prespecified and not fitted; failure to outperform is informative and should not be hidden.

## Reproduction

See `command.txt`, `manifest.json`, `results.json`, and the cohort CSV files in this directory.
