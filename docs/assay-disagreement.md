# Assay disagreement prediction protocol

## Scientific question

Can structure and evidence available before a new experiment predict whether the independent later study will disagree with prior exact-IC50 evidence? Disagreement is a reproducibility outcome, not evidence of drug failure, efficacy, toxicity, or clinical relevance.

## Cohort and outcome

The five prespecified human single-protein targets and normalized ChEMBL inputs from the multi-target benchmark are reused. Within each target, only the largest compatible target/organism/IC50/assay-type/BAO context is retained. Relations must be exact and units must have passed normalization. Records carrying a ChEMBL data-validity comment, a potential-duplicate flag, or a computed/reported pChEMBL mismatch greater than the pipeline tolerance are explicitly excluded and counted.

For a molecule with at least two publication years, all measurements in its latest year form the held-out event. Earlier measurements form its history. Activity IDs, ChEMBL document IDs, and assay IDs must be disjoint. The primary binary outcome is an absolute difference of at least 1.0 log10 unit between the earlier-evidence median and latest-year median. This ten-fold threshold was fixed before evaluation. Crossing pActivity 6.0 is a descriptive secondary outcome.

## Leakage controls

- Development outcomes have held-out publication years through 2019; evaluation outcomes start in 2020.
- A molecule appearing in evaluation for any target is removed from development for every target.
- Predictors contain only RDKit descriptors, target identity, and evidence/provenance summaries computed from strictly earlier records.
- Held-out assay descriptions, identifiers, documents, values, and outcome-derived fields are prohibited predictors.
- NPASS/COCONUT status is used only for subgroup evaluation, not prediction, because current database releases postdate some held-out studies.
- Model definitions and regularization are fixed. The untouched evaluation set is not used for tuning or selection.
- Novel-scaffold evaluation is reported as a sensitivity analysis.

## Models and uncertainty

The primary model is an interpretable L2 logistic regression using molecular descriptors, earlier-study consistency, and provenance counts. Fixed comparators include development prevalence, target only, molecular only, prior potency, history only, provenance only, and group-removal ablations.

Metrics are AUROC, PR-AUC, Brier score, log loss, ten-bin expected calibration error, calibration intercept, and calibration slope. Confidence intervals and paired model differences resample held-out publication clusters within target. Per-compound risk intervals refit the primary model on target-stratified development bootstrap samples. These intervals reflect development-sample instability, not all sources of uncertainty.

The primary model is called supported only if its study-clustered AUROC 95% interval lies above 0.5 and its PR-AUC interval lies above the untouched evaluation prevalence. Failure is reported directly; a serialized failed model remains an audit artifact, not a production predictor.

## Reproduction

```powershell
cd backend
python -m app.pipeline assay-disagreement `
  --specs ../experiments/multitarget_targets.json `
  --output-dir ../experiments/assay_disagreement `
  --bootstrap-iterations 1000 `
  --prediction-interval-iterations 200 `
  --seed 2026
```

Outputs include the provenance-preserving event dataset, predictions, clustered metrics and comparisons, coefficients, held-out study clusters, JSON results, input hashes, a plot, a generated report, and an explicitly labeled research model artifact.

## Limitations

ChEMBL does not encode substrate concentration, incubation duration, construct, detection method, and laboratory consistently enough to model detailed protocol transfer. Publication is only a proxy for independent study. A large multi-compound publication can dominate event counts, which is why inference is study-clustered. Five targets are insufficient for a broad cross-target claim, cross-organism transfer is not estimable, and the natural-product subgroup is expected to be underpowered.
