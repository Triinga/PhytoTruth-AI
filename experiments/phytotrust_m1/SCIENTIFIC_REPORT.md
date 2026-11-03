# PhytoTRUST-M1 development-only scientific evaluation

## Scope

M1 used only the CHEMBL222 development artifact (966 measurements, 877 compounds, 61 publications) and full ChEMBL metadata for 961 referenced assay instances. The future-evaluation history, index, and quarantined labels were not opened.

The complete frozen target snapshot had already been exhausted. No post-outcome compatibility relaxation or outcome-guided measurement acquisition was performed. The targeted data improvement was development-only assay-metadata retrieval and literal protocol extraction.

## Method

M1 predicts a new-study measurement as the historical compound median plus a regularized evidence-conditioned correction and a publication random intercept. Its Student-t residual scale is conditioned on historical dispersion, support, protocol missingness/diversity, and history span. The unseen publication effect is integrated out. Held-out protocol fields are never predictors.

The individual ingredients are established statistical ideas. The tested contribution is their timing-safe integration for new-study evidence transferability with explicit evidence-conditioned uncertainty; this experiment does not establish novelty.

## Protocol metadata coverage

Mean direct-field completeness was 0.497; mean missing fraction across the prespecified protocol fields was 0.909. Literal protocol categories were found for 907/907 compound-publication edges. Unknown values remained unknown.

## Results

| Validation | Scope | Method | n | MAE | RMSE | NLL | 80% coverage / width | 90% coverage / width | calibration error |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| forward_temporal | all | m1_homoscedastic | 25 | 0.519 | 0.686 | 1.510 | 0.760 / 1.196 | 0.840 / 1.448 | 0.050 |
| forward_temporal | natural_products | m1_homoscedastic | 3 | 0.870 | 1.222 | 2.424 | 0.667 / 1.173 | 0.667 / 1.551 | 0.183 |
| forward_temporal | all | m1_no_protocol | 25 | 0.518 | 0.685 | 1.506 | 0.760 / 1.196 | 0.840 / 1.448 | 0.050 |
| forward_temporal | natural_products | m1_no_protocol | 3 | 0.870 | 1.222 | 2.424 | 0.667 / 1.173 | 0.667 / 1.551 | 0.183 |
| forward_temporal | all | m1_no_study_variance | 25 | 0.515 | 0.667 | 1.364 | 0.640 / 5.911 | 0.760 / 8.218 | 0.150 |
| forward_temporal | natural_products | m1_no_study_variance | 3 | 0.847 | 1.149 | 2.085 | 0.667 / 1.533 | 0.667 / 2.132 | 0.183 |
| forward_temporal | all | median_potency | 25 | 0.505 | 0.669 | 1.167 | 0.560 / 0.879 | 0.680 / 1.222 | 0.230 |
| forward_temporal | natural_products | median_potency | 3 | 0.854 | 1.176 | 2.263 | 0.667 / 1.208 | 0.667 / 1.680 | 0.183 |
| forward_temporal | all | phytotrust_m0 | 25 | 0.526 | 0.701 | 1.145 | 0.680 / 1.149 | 0.840 / 1.455 | 0.090 |
| forward_temporal | natural_products | phytotrust_m0 | 3 | 0.895 | 1.244 | 2.468 | 0.667 / 1.174 | 0.667 / 1.554 | 0.183 |
| forward_temporal | all | phytotrust_m1 | 25 | 0.539 | 0.724 | 1.521 | 0.600 / 5.888 | 0.760 / 8.157 | 0.170 |
| forward_temporal | natural_products | phytotrust_m1 | 3 | 0.919 | 1.276 | 2.480 | 0.667 / 1.312 | 0.667 / 1.761 | 0.183 |
| forward_temporal | all | student_t_regression | 25 | 0.517 | 0.690 | 1.100 | 0.640 / 1.027 | 0.840 / 1.428 | 0.110 |
| forward_temporal | natural_products | student_t_regression | 3 | 0.880 | 1.220 | 2.109 | 0.667 / 1.142 | 0.667 / 1.588 | 0.183 |
| forward_temporal | all | xgboost | 25 | 0.504 | 0.668 | 1.054 | 0.640 / 1.050 | 0.720 / 1.461 | 0.170 |
| forward_temporal | natural_products | xgboost | 3 | 0.854 | 1.182 | 2.165 | 0.667 / 1.474 | 0.667 / 2.050 | 0.183 |
| study_grouped | all | m1_homoscedastic | 59 | 0.460 | 0.632 | 0.990 | 0.729 / 1.226 | 0.847 / 1.649 | 0.062 |
| study_grouped | natural_products | m1_homoscedastic | 6 | 0.848 | 1.170 | 2.150 | 0.667 / 1.242 | 0.667 / 1.667 | 0.183 |
| study_grouped | all | m1_no_protocol | 59 | 0.455 | 0.619 | 0.970 | 0.763 / 1.273 | 0.847 / 1.724 | 0.045 |
| study_grouped | natural_products | m1_no_protocol | 6 | 0.846 | 1.150 | 2.194 | 0.667 / 1.293 | 0.667 / 1.749 | 0.183 |
| study_grouped | all | m1_no_study_variance | 59 | 0.455 | 0.623 | 1.166 | 0.678 / 1.971 | 0.797 / 2.741 | 0.113 |
| study_grouped | natural_products | m1_no_study_variance | 6 | 0.808 | 1.123 | 2.341 | 0.667 / 2.101 | 0.667 / 2.922 | 0.183 |
| study_grouped | all | median_potency | 59 | 0.452 | 0.623 | 1.033 | 0.797 / 1.425 | 0.864 / 1.981 | 0.019 |
| study_grouped | natural_products | median_potency | 6 | 0.854 | 1.176 | 2.159 | 0.667 / 1.266 | 0.667 / 1.761 | 0.183 |
| study_grouped | all | phytotrust_m0 | 59 | 0.453 | 0.627 | 0.973 | 0.780 / 1.364 | 0.898 / 1.805 | 0.011 |
| study_grouped | natural_products | phytotrust_m0 | 6 | 0.856 | 1.179 | 2.197 | 0.667 / 1.374 | 0.667 / 1.812 | 0.183 |
| study_grouped | all | phytotrust_m1 | 59 | 0.511 | 0.699 | 1.372 | 0.695 / 4.395 | 0.797 / 6.052 | 0.104 |
| study_grouped | natural_products | phytotrust_m1 | 6 | 0.882 | 1.254 | 2.861 | 0.667 / 2.163 | 0.667 / 2.946 | 0.183 |
| study_grouped | all | student_t_regression | 59 | 0.453 | 0.626 | 0.964 | 0.780 / 1.333 | 0.881 / 1.854 | 0.019 |
| study_grouped | natural_products | student_t_regression | 6 | 0.857 | 1.181 | 1.979 | 0.667 / 1.311 | 0.667 / 1.823 | 0.183 |
| study_grouped | all | xgboost | 59 | 0.450 | 0.621 | 1.049 | 0.831 / 1.748 | 0.915 / 2.431 | 0.023 |
| study_grouped | natural_products | xgboost | 6 | 0.851 | 1.171 | 2.032 | 0.667 / 1.522 | 0.667 / 2.117 | 0.183 |

## Primary interpretation

Study-grouped M1 NLL was 1.372; the strongest comparator was student_t_regression at 0.964. Forward-temporal M1 NLL was 1.521; the strongest comparator was xgboost at 1.054.

M1-minus-student_t_regression study-grouped NLL 95% CI: [0.0730325031618373, 0.8964485528806695]. M1-minus-xgboost forward-temporal NLL 95% CI: [0.016751565270393835, 1.6300405157611504]. Negative values favor M1.

## Ablations

Study-grouped NLL: full M1 1.372; without historical protocol features 0.970; homoscedastic M1 0.990; without study variance 1.166; M0 0.973.

Forward-temporal NLL: full M1 1.521; without protocol 1.506; homoscedastic 1.510; without study variance 1.364; M0 1.145.

An apparent gain is attributable to the proposed protocol-conditioned scale only if full M1 improves over `m1_no_protocol` with uncertainty excluding zero in both validation schemes. The same logic applies to heteroscedasticity and study variance.

## Natural products and connectivity

Only 6 connected natural-product transfer events were available in the full development graph. NP results are exploratory. The multiple-prior-study stratum contains only 3 study-grouped and 1 temporal events, so improved performance cannot be attributed reliably to deeper compound connectivity.

## Limitations

- M1 predicts measurement transferability, not true bioactivity, efficacy, safety, or clinical utility.
- Publication is an imperfect proxy for independent study and laboratory.
- ChEMBL protocol descriptions are retrospective records. Only historical descriptions were used, but their original pre-experiment availability is not proven.
- Protocol fields are sparse and literal text extraction cannot recover unreported conditions.
- Leave-publication-out compound residuals are crossed and dependent; publication bootstrap intervals remain approximate.
- The study-grouped analysis is retrospective; only forward-temporal validation obeys a strict time filtration.
- One target, 29 repeated compounds, and very small NP strata prevent broad generalization.

## Publication-readiness verdict

M1 did not meet the frozen development-only support rule. The evidence-conditioned formulation remains a valid negative feasibility result, but the available data do not support a claim of superiority, methodological novelty, or IEEE journal readiness. M1 must not be promoted to the quarantined final test.
