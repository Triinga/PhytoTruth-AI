# PhytoTRUST-M0 on the study-connected CHEMBL222 dataset

## Scope and leakage contract

This experiment used only `development_measurements.csv`: 966 real development measurements, 877 compounds, and 61 ChEMBL publications. No future-evaluation history or quarantined label file was opened.

Study-grouped validation leaves out entire publications and excludes their activity and assay identifiers. It is a retrospective domain-holdout analysis and is not described as prospective because other development publications can be later in time. Forward-temporal validation uses only strictly earlier publication years. Same-compound historical measurements are intentionally allowed because prediction for a previously measured compound is the declared estimand; the held-out publication's records are never used in its features or fitting data.

M0 uses the prespecified Student-t residual model with a normally distributed publication random intercept. The random effect for a new publication is integrated out. The non-hierarchical ablation has the same location and Student-t likelihood but fixes between-study variance to zero.

## Validation metrics

| Validation | Scope | Method | n | studies | MAE | RMSE | NLL | 80% coverage / width | 90% coverage / width |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| forward_temporal | all | global_mean | 25 | 8 | 0.865 | 1.018 | 1.526 | 0.800 / 2.566 | 1.000 / 3.568 |
| forward_temporal | natural_products | global_mean | 3 | 3 | 1.046 | 1.216 | 1.767 | 0.667 / 2.622 | 1.000 / 3.646 |
| forward_temporal | one_prior_study | global_mean | 24 | 7 | 0.843 | 1.000 | 1.503 | 0.833 / 2.564 | 1.000 / 3.565 |
| forward_temporal | multiple_prior_studies | global_mean | 1 | 1 | 1.395 | 1.395 | 2.100 | 0.000 / 2.620 | 1.000 / 3.643 |
| forward_temporal | all | median_potency | 25 | 8 | 0.505 | 0.669 | 1.167 | 0.560 / 0.879 | 0.680 / 1.222 |
| forward_temporal | natural_products | median_potency | 3 | 3 | 0.854 | 1.176 | 2.263 | 0.667 / 1.208 | 0.667 / 1.680 |
| forward_temporal | one_prior_study | median_potency | 24 | 7 | 0.503 | 0.673 | 1.174 | 0.542 / 0.860 | 0.667 / 1.195 |
| forward_temporal | multiple_prior_studies | median_potency | 1 | 1 | 0.554 | 0.554 | 0.997 | 1.000 / 1.340 | 1.000 / 1.863 |
| forward_temporal | all | nonhierarchical_student_t | 25 | 8 | 0.517 | 0.690 | 1.100 | 0.640 / 1.027 | 0.840 / 1.428 |
| forward_temporal | natural_products | nonhierarchical_student_t | 3 | 3 | 0.880 | 1.220 | 2.109 | 0.667 / 1.142 | 0.667 / 1.588 |
| forward_temporal | one_prior_study | nonhierarchical_student_t | 24 | 7 | 0.516 | 0.695 | 1.102 | 0.625 / 1.021 | 0.833 / 1.420 |
| forward_temporal | multiple_prior_studies | nonhierarchical_student_t | 1 | 1 | 0.547 | 0.547 | 1.058 | 1.000 / 1.164 | 1.000 / 1.618 |
| forward_temporal | all | phytotrust_m0_connected | 25 | 8 | 0.526 | 0.701 | 1.145 | 0.680 / 1.149 | 0.840 / 1.455 |
| forward_temporal | natural_products | phytotrust_m0_connected | 3 | 3 | 0.895 | 1.244 | 2.468 | 0.667 / 1.174 | 0.667 / 1.554 |
| forward_temporal | one_prior_study | phytotrust_m0_connected | 24 | 7 | 0.525 | 0.707 | 1.153 | 0.667 / 1.149 | 0.833 / 1.452 |
| forward_temporal | multiple_prior_studies | phytotrust_m0_connected | 1 | 1 | 0.545 | 0.545 | 0.946 | 1.000 / 1.152 | 1.000 / 1.534 |
| study_grouped | all | global_mean | 59 | 18 | 0.911 | 1.057 | 1.548 | 0.797 / 2.883 | 1.000 / 4.009 |
| study_grouped | natural_products | global_mean | 6 | 6 | 1.417 | 1.535 | 2.125 | 0.333 / 2.891 | 1.000 / 4.020 |
| study_grouped | one_prior_study | global_mean | 56 | 17 | 0.909 | 1.059 | 1.548 | 0.786 / 2.883 | 1.000 / 4.009 |
| study_grouped | multiple_prior_studies | global_mean | 3 | 3 | 0.953 | 1.021 | 1.540 | 1.000 / 2.885 | 1.000 / 4.012 |
| study_grouped | all | median_potency | 59 | 18 | 0.452 | 0.623 | 1.033 | 0.797 / 1.425 | 0.864 / 1.981 |
| study_grouped | natural_products | median_potency | 6 | 6 | 0.854 | 1.176 | 2.159 | 0.667 / 1.266 | 0.667 / 1.761 |
| study_grouped | one_prior_study | median_potency | 56 | 17 | 0.449 | 0.627 | 1.045 | 0.786 / 1.428 | 0.857 / 1.985 |
| study_grouped | multiple_prior_studies | median_potency | 3 | 3 | 0.495 | 0.546 | 0.801 | 1.000 / 1.369 | 1.000 / 1.904 |
| study_grouped | all | nonhierarchical_student_t | 59 | 18 | 0.453 | 0.626 | 0.964 | 0.780 / 1.333 | 0.881 / 1.854 |
| study_grouped | natural_products | nonhierarchical_student_t | 6 | 6 | 0.857 | 1.181 | 1.979 | 0.667 / 1.311 | 0.667 / 1.823 |
| study_grouped | one_prior_study | nonhierarchical_student_t | 56 | 17 | 0.450 | 0.630 | 0.966 | 0.786 / 1.335 | 0.875 / 1.856 |
| study_grouped | multiple_prior_studies | nonhierarchical_student_t | 3 | 3 | 0.497 | 0.549 | 0.932 | 0.667 / 1.300 | 1.000 / 1.808 |
| study_grouped | all | phytotrust_m0_connected | 59 | 18 | 0.453 | 0.627 | 0.973 | 0.780 / 1.364 | 0.898 / 1.805 |
| study_grouped | natural_products | phytotrust_m0_connected | 6 | 6 | 0.856 | 1.179 | 2.197 | 0.667 / 1.374 | 0.667 / 1.812 |
| study_grouped | one_prior_study | phytotrust_m0_connected | 56 | 17 | 0.450 | 0.630 | 0.977 | 0.786 / 1.364 | 0.893 / 1.806 |
| study_grouped | multiple_prior_studies | phytotrust_m0_connected | 3 | 3 | 0.506 | 0.554 | 0.910 | 0.667 / 1.356 | 1.000 / 1.792 |

## Primary comparison

Study-grouped M0: MAE 0.453, RMSE 0.627, NLL 0.973. Median potency: MAE 0.452, RMSE 0.623, NLL 1.033. Matched no-study-variance ablation NLL: 0.964.

Forward-temporal M0: MAE 0.526, RMSE 0.701, NLL 1.145. Median potency: MAE 0.505, RMSE 0.669, NLL 1.167.

The study-cluster bootstrap 95% interval for study-grouped M0-minus-median NLL is [-0.11521287361477824, 0.03558947719903971]. Negative values favor M0. These intervals use publication as the resampling unit and are not independent confirmation.

## Ablations and dominant-publication sensitivity

The matched non-hierarchical Student-t ablation removes only the publication random-effect variance. Median potency removes the learned global correction and uses a single empirical residual scale. Global mean removes compound history as well. M0 was not significantly better than median potency or the matched no-study-variance ablation in the publication-cluster bootstrap.

The largest connected publication was CHEMBL1156106 (12 held-out transfer events). After excluding those held-out predictions, study-grouped M0 NLL was 1.083, versus 1.065 for the no-study-variance ablation. Refitting the final variance model without that publication changed between-study SD from 0.3235 to 0.2889 pChEMBL.

## Connectivity analysis

The variance-estimation subgraph contains 59 edges from 29 repeated compounds and 18 publications. Only 9 publications have at least two shared-compound anchors, 6 have at least three, and 1 compound appears in at least three studies.

Metrics stratified by one versus multiple prior studies and by single- versus multi-bridge publications are included in `metrics.csv`. The multiple-prior stratum is too small for a stable claim; differences are descriptive rather than causal evidence that connectivity improves prediction.

## Parameter uncertainty and identifiability

Between-study variance is not reliably identifiable: the connected graph and/or fold/bootstrap stability criteria failed. The estimate is a regularized feasibility quantity, not a definitive variance component.

Final between-study SD: 0.3235 pChEMBL; cross-fold range: [0.23581117493312012, 0.41672749333637854]. Model-based parametric-bootstrap 95% interval: [0.0743617719419368, 0.6004446663918986]. The bootstrap is conditional on M0 and cannot account for publication selection, dependence induced by leave-publication-out medians, or database curation uncertainty.

## Natural products

Only 10 verified natural-product compounds occur in the full development input, and 6 appear in at least one validation event. NP status is used only for subgroup reporting, never as a predictor or fitting weight. Subgroup estimates are exploratory.

## Limitations

- M0 predicts a comparable experimental measurement, not true potency, efficacy, safety, or clinical benefit.
- Publication is a proxy for study/laboratory and may merge or split actual experimental domains.
- Leave-publication-out residual events share compound measurements and are statistically dependent; the study-cluster bootstrap does not remove all crossed dependence.
- Study-grouped validation is not temporally prospective; only the forward-temporal analysis enforces the time filtration.
- Sparse repeated compounds, a single target/context, and very few natural products limit generalization.
- A structural graph gate is necessary but not sufficient for variance identifiability.

## Scientific verdict

Development-only evidence does not establish M0 as a superior or statistically identified foundation. The connected graph permits a more relevant feasibility analysis, but sparse anchors and/or inconsistent held-out performance require more connected independent studies before M1 or a journal-level claim.
