# Evidence Reliability Engine scientific investigation

Generated 2026-10-20T14:15:03.929154+00:00 from normalized ChEMBL measurements and the recorded NPASS structure cross-match.

## Prespecified validation design

The latest-publication-per-molecule protocol was retained: all score inputs precede each molecule's held-out publication evidence, and documents are disjoint across time. Molecules were then deterministically partitioned, using only a hash of molecule ID and NPASS membership, into 98 development and 64 untouched-for-revision evaluation molecules. The evaluation outcomes were not inputs to the selection gate. This partition is untouched for this revision procedure, but it is not an external validation dataset because aggregate results on the parent cohort were reported previously.

## Why the current heuristic fails

The current score subtracts inconsistency once in `evidence_score`, then adds a composite reliability penalty that includes the same MAD/contradiction signal, and finally subtracts bootstrap interval width. This double/triple use of dispersion is not calibrated to future assay outcomes. Its positive support term counts assays, while independent publications are only weakly and saturatingly represented inside the composite reliability value. Consequently, the heuristic can reverse a robust median ordering when sparse or heterogeneous evidence produces a large penalty.

On the untouched evaluation partition, current PhytoTruth minus median potency Kendall tau-b was -0.057 (95% CI -0.151 to 0.006; p=0.1639). nDCG@10 difference was -0.044 (95% CI -0.089 to 0.104).

Pairwise error analysis found 72 harmful and 14 beneficial changes to median-correct orderings among 2007 comparable evaluation pairs. Dominant harmful components were assay_support: 7, experimental_consistency: 65. See `compound_error_analysis.csv` and `harmful_reversals.csv` for traceable molecule-level details.

## Natural products

The untouched NPASS-verified evaluation subgroup has only 12 molecules. Current PhytoTruth minus median nDCG@10 was -0.409 (95% CI -0.685 to 0.000); Kendall difference was -0.212. Median nDCG was 0.850, versus 0.441 for the current heuristic. This sample is too small for a reliable subgroup superiority claim. COCONUT did not contribute verified compounds to this target artifact; this subgroup is NPASS-only.

## Development-only decision and defensible change

The prespecified gate was: Select the current evidence adjustment only when development-set nDCG@K and Kendall tau-b both improve over median potency, both paired 95% interval lower bounds exceed zero, and both two-sided score-swap p-values are below 0.05. It did not pass on development data, so the selected ranking is `median_potency`. No new weights were fitted. The scientifically defensible change is therefore a conservative guardrail: use robust median potency for prioritization and retain reliability fields as descriptive warnings, not as validated ranking probabilities or penalties.

Untouched evaluation median metrics were nDCG@10=0.915, precision@10=1.000, enrichment=1.730, and Kendall tau-b=0.887. Current-heuristic metrics were nDCG@10=0.871, precision@10=1.000, enrichment=1.730, and Kendall tau-b=0.830.

## Predictive-value and ablation interpretation

`component_predictive_value.csv` reports development and evaluation Spearman associations, molecule-bootstrap 95% intervals, and permutation p-values for two questions: whether a component points in the direction of later potency drift, and whether its risk signal predicts absolute future error. These are association tests, not probability calibration. `metric_summary.csv` compares maximum, mean, most-recent, lower-confidence-bound, and median potency with every existing ablation.

## Limitations and required data

- The compatible ChEMBL context is one target/context; generalization to other targets is unknown.
- The revision evaluation partition contains 64 molecules and only 12 NPASS-verified natural products, so natural-product confidence intervals are wide.
- Publication identifiers are a proxy for independent studies; laboratory, protocol, batch, and replicate identifiers are unavailable.
- Bootstrap intervals over repeated measurements do not capture systematic assay bias.
- A scientifically credible learned reliability adjustment needs a larger multi-target dataset, externally held-out publications/laboratories, richer assay-protocol covariates, and substantially more COCONUT/NPASS-verified molecules. Until then, reliability indicators should remain explanatory warnings.
