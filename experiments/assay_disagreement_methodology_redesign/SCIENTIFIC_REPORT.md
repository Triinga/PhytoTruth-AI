# Assay Disagreement Methodology Redesign

## Scope and leakage controls

This development-only analysis used 212 pre-2020 events (50 disagreements) from 141 publication groups and 5 targets. It excluded 319 rows in the locked 2020+ partition. The separately assembled 90-case and 94-case independent datasets were not loaded.

The primary outcome is absolute held-out pChEMBL shift. Binary ten-fold disagreement is secondary. Every validation fold separates molecules; study-grouped folds also separate publication groups, forward folds train strictly before their test years, and leave-target-out folds hold out a complete target.

## Why the previous method failed

The original development set has 112 singleton publication groups out of 141; consequently, study-specific shifts cannot be learned robustly. The earlier investigation also measured study-level label clustering (ICC 0.288) and found that its calibrated candidates were worse than a constant-prevalence model in development. Protocol enrichment did not recover the central missing conditions: across 1,050 assays, cell type and strain were never structured, while temperature, pH, methodology, endpoint, and parameter types were reported or detected in only 2.8%, 4.2%, 12.6%, 11.7%, and 1.2% respectively. Reporting indicators therefore cannot substitute for actual protocol values.

## Development validation

| Scheme | Scope | Method | n | AUROC | PR-AUC | MAE | Brier |
|---|---|---|---:|---:|---:|---:|---:|
| study_grouped | all | constant | 212 | 0.408 | 0.203 | 0.527 | 0.183 |
| study_grouped | all | historical_beta | 212 | 0.558 | 0.272 | 0.732 | 0.180 |
| study_grouped | all | legacy_logistic | 212 | 0.541 | 0.262 | — | 0.211 |
| study_grouped | all | huber_magnitude | 212 | 0.567 | 0.289 | 0.595 | — |
| study_grouped | all | rf_history | 212 | 0.575 | 0.266 | 0.593 | 0.204 |
| study_grouped | all | rf_protocol | 212 | 0.572 | 0.283 | 0.586 | 0.196 |
| study_grouped | natural_products | constant | 47 | 0.396 | 0.261 | 0.512 | 0.215 |
| study_grouped | natural_products | historical_beta | 47 | 0.549 | 0.376 | 0.788 | 0.212 |
| study_grouped | natural_products | legacy_logistic | 47 | 0.548 | 0.359 | — | 0.250 |
| study_grouped | natural_products | huber_magnitude | 47 | 0.569 | 0.378 | 0.658 | — |
| study_grouped | natural_products | rf_history | 47 | 0.682 | 0.408 | 0.560 | 0.213 |
| study_grouped | natural_products | rf_protocol | 47 | 0.547 | 0.382 | 0.528 | 0.219 |
| forward_temporal | all | constant | 131 | 0.552 | 0.235 | 0.487 | 0.171 |
| forward_temporal | all | historical_beta | 131 | 0.544 | 0.243 | 0.710 | 0.172 |
| forward_temporal | all | legacy_logistic | 131 | 0.500 | 0.224 | — | 0.226 |
| forward_temporal | all | huber_magnitude | 131 | 0.524 | 0.271 | 0.801 | — |
| forward_temporal | all | rf_history | 131 | 0.536 | 0.275 | 0.628 | 0.231 |
| forward_temporal | all | rf_protocol | 131 | 0.574 | 0.346 | 0.572 | 0.200 |
| forward_temporal | natural_products | constant | 35 | 0.554 | 0.362 | 0.576 | 0.231 |
| forward_temporal | natural_products | historical_beta | 35 | 0.563 | 0.442 | 0.923 | 0.235 |
| forward_temporal | natural_products | legacy_logistic | 35 | 0.543 | 0.366 | — | 0.292 |
| forward_temporal | natural_products | huber_magnitude | 35 | 0.496 | 0.380 | 0.826 | — |
| forward_temporal | natural_products | rf_history | 35 | 0.582 | 0.403 | 0.642 | 0.254 |
| forward_temporal | natural_products | rf_protocol | 35 | 0.531 | 0.430 | 0.656 | 0.251 |
| leave_target_out | all | constant | 212 | 0.413 | 0.205 | 0.528 | 0.182 |
| leave_target_out | all | historical_beta | 212 | 0.558 | 0.272 | 0.732 | 0.180 |
| leave_target_out | all | legacy_logistic | 212 | 0.539 | 0.264 | — | 0.225 |
| leave_target_out | all | huber_magnitude | 212 | 0.487 | 0.238 | 0.670 | — |
| leave_target_out | all | rf_history | 212 | 0.579 | 0.286 | 0.561 | 0.200 |
| leave_target_out | all | rf_protocol | 212 | 0.515 | 0.263 | 0.562 | 0.202 |
| leave_target_out | natural_products | constant | 47 | 0.300 | 0.229 | 0.523 | 0.216 |
| leave_target_out | natural_products | historical_beta | 47 | 0.549 | 0.376 | 0.788 | 0.212 |
| leave_target_out | natural_products | legacy_logistic | 47 | 0.416 | 0.273 | — | 0.306 |
| leave_target_out | natural_products | huber_magnitude | 47 | 0.368 | 0.255 | 0.786 | — |
| leave_target_out | natural_products | rf_history | 47 | 0.600 | 0.360 | 0.581 | 0.242 |
| leave_target_out | natural_products | rf_protocol | 47 | 0.544 | 0.353 | 0.525 | 0.225 |

## Scientific interpretation

Magnitude modeling was evaluated because it retains more information than binary thresholding. In study-grouped validation, the protocol RF had AUROC 0.572, PR-AUC 0.283, and MAE 0.586; the constant baseline MAE was 0.527. For natural products, protocol-RF AUROC was 0.5465367965367965 versus 0.6818181818181819 for history-only RF. These are development estimates, not independent validation or evidence of deployment readiness.

In the paired study-cluster bootstrap, adding molecular and prior protocol-completeness features changed overall AUROC by -0.002 (95% CI [-0.10287934069188079, 0.09290764360018088]) and MAE by -0.007 (95% CI [-0.04531318438869049, 0.03482457541735367]; negative favors the richer model). For natural products, the AUROC change was -0.135 (95% CI [-0.3086449430199431, 0.02326035610465109]). These intervals determine whether the added feature group has defensible value; point estimates alone do not.

**Readiness decision: `not_ready_for_new_independent_test`.** At least one prespecified generalization criterion failed; additional independent test use is not scientifically justified yet.

Random-forest tree quantiles are descriptive ensemble intervals, not validated confidence or prediction intervals. The historical-beta score is a Jeffreys-smoothed prior contradiction frequency, not a validated probability. Natural-product status was never used as a model feature.

At the fixed 0.5 descriptive RF vote threshold, the study-grouped natural-product analysis missed 12 disagreements and produced 1 false alarms. Compound-level provenance and residuals are in `natural_product_error_analysis.csv`; this threshold was not tuned.

## Limitations

- The development set is small, with only 50 disagreements and 14 natural-product disagreements; subgroup confidence intervals are consequently wide.
- Protocol features mainly encode whether metadata were reported, not the missing experimental conditions themselves.
- Leave-target-out validation is severe with only five targets, while forward-temporal folds cover only the later part of the development era.
- Disagreement indicates a large difference between comparable database measurements; it does not establish which experiment is correct or imply clinical failure.
- No result here updates or reinterprets any previous independent evaluation.
