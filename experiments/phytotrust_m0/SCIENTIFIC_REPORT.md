# PhytoTRUST-M0 scientific validation

## Scope

This experiment used only the 212 real pre-2020 development events from 141 ChEMBL publication groups and 5 targets. The locked 2020+ partition and previously examined independent cohorts were not used for fitting, selection, or evaluation.

M0 models held-out pChEMBL as the historical median plus a target correction, an independent-study random intercept, and Student-t residual noise. For blind transfer, the held-out study effect is integrated over its learned target-specific population distribution. No outcome from the held-out publication is used at prediction time.

## Results

| Validation | Scope | Method | n | MAE | RMSE | NLL | 80% coverage / width | 95% coverage / width |
|---|---|---|---:|---:|---:|---:|---:|---:|
| forward_temporal | all | global_mean | 131 | 1.084 | 1.292 | 1.729 | 0.878 / 3.620 | 0.992 / 6.555 |
| forward_temporal | natural_products | global_mean | 35 | 1.055 | 1.269 | 1.715 | 0.829 / 3.615 | 1.000 / 6.546 |
| forward_temporal | all | median_potency | 131 | 0.617 | 0.956 | 1.262 | 0.824 / 2.288 | 0.947 / 4.143 |
| forward_temporal | natural_products | median_potency | 35 | 0.757 | 1.086 | 1.449 | 0.743 / 2.287 | 0.914 / 4.141 |
| forward_temporal | all | nonhierarchical_student_t | 131 | 0.676 | 0.992 | 1.312 | 0.802 / 2.273 | 0.954 / 4.115 |
| forward_temporal | natural_products | nonhierarchical_student_t | 35 | 0.749 | 1.099 | 1.442 | 0.743 / 2.274 | 0.943 / 4.119 |
| forward_temporal | all | phytotrust_m0 | 131 | 0.710 | 1.022 | 1.410 | 0.817 / 2.402 | 0.931 / 3.897 |
| forward_temporal | natural_products | phytotrust_m0 | 35 | 0.749 | 1.108 | 1.468 | 0.743 / 2.077 | 0.886 / 3.459 |
| study_grouped | all | global_mean | 212 | 1.115 | 1.360 | 1.776 | 0.783 / 3.566 | 0.995 / 6.457 |
| study_grouped | natural_products | global_mean | 47 | 1.092 | 1.315 | 1.748 | 0.766 / 3.562 | 1.000 / 6.449 |
| study_grouped | all | median_potency | 212 | 0.681 | 1.058 | 1.353 | 0.764 / 2.037 | 0.910 / 3.688 |
| study_grouped | natural_products | median_potency | 47 | 0.684 | 0.989 | 1.351 | 0.702 / 2.021 | 0.915 / 3.659 |
| study_grouped | all | nonhierarchical_student_t | 212 | 0.716 | 1.074 | 1.378 | 0.797 / 2.218 | 0.920 / 4.017 |
| study_grouped | natural_products | nonhierarchical_student_t | 47 | 0.672 | 0.973 | 1.330 | 0.787 / 2.221 | 0.915 / 4.022 |
| study_grouped | all | phytotrust_m0 | 212 | 0.726 | 1.081 | 1.429 | 0.821 / 2.435 | 0.915 / 4.006 |
| study_grouped | natural_products | phytotrust_m0 | 47 | 0.672 | 0.978 | 1.289 | 0.787 / 2.084 | 0.915 / 3.542 |

## Interpretation

In study-grouped validation, M0 had MAE 0.726, RMSE 1.081, and NLL 1.429. The matched non-hierarchical Student-t model had MAE 0.716, RMSE 1.074, and NLL 1.378. Equal point errors are expected because blind integration of a zero-mean unseen-study effect changes uncertainty, not the predictive mean.

The cluster-bootstrap comparison intervals are stored in `bootstrap_differences.csv` and `results.json`. A negative M0-minus-baseline interval favors M0. Results are development-only and do not establish novelty or independent generalization.

M0 did not improve the primary probabilistic result: its study-grouped NLL was 1.429, versus 1.353 for median potency and 1.378 for the matched non-hierarchical model. In forward-time validation, the corresponding values were 1.410, 1.262, and 1.312. The natural-product study-grouped point estimate favored M0 on NLL, but its paired study-bootstrap intervals included no effect and the direction did not reproduce temporally.

Removing the largest publication (CHEMBL5725214, 25 events) changed study-grouped M0 NLL to 1.461; the non-hierarchical value became 1.443. The main negative conclusion therefore does not depend solely on that publication.

## Variance identifiability

Target-specific between-study variance is not statistically identifiable with the available cluster structure. M0 estimates are regularized feasibility estimates, not definitive variance components.

Targets meeting the structural identifiability gate: 0/5. Partial pooling stabilizes estimation but cannot create information absent from singleton studies.

## Natural products

The verified-natural-product subgroup contains 47 events. It was evaluation-only: NP status was not a predictor or fitting weight. Its estimates are exploratory because study and target coverage are sparse.

## Limitations

- M0 estimates a measurement distribution, not true potency, efficacy, safety, or which study is correct.
- The dataset contains event-level study medians rather than every within-assay replicate, limiting separation of within-assay and between-study noise.
- Target-specific variance is weakly identified because 112/141 studies are singletons and only 10 studies have at least three events.
- Student-t degrees of freedom and pooling scales were prespecified from the methodology, not tuned on held-out outcomes.
- Temporal validation is compound-disjoint and therefore tests a strict form of transfer; it may omit historically repeated compounds from training by design.
- This valid M0 implementation is a feasibility baseline, not evidence that the method is novel or publication-ready.

## Scientific assessment

M0 is a valid, reproducible hierarchical feasibility baseline, but the available development evidence does not support it as a journal-level methodological foundation yet. It failed to improve NLL over median potency or the matched non-hierarchical model, temporal performance was worse, and target-specific study variance did not pass the identifiability gate. More multi-compound independent studies are required before extending the hierarchy.
