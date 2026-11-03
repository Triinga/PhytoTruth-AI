# Expanded independent validation of frozen disagreement models

## Locked protocol

The 94-event assay-disagreement expansion was evaluated only after the target,
measurement, natural-product, and disagreement-label construction had finished.
The validation protocol was recorded before frozen-model scoring. Neither model
was retrained, recalibrated, tuned, or selected using these outcomes.

The evaluated artifacts are the original interpretable full logistic model and
the study-generalized artifact whose previously selected fallback is the
constant development prevalence (`0.2358490566`). SHA-256 checksums were
identical before and after inference:

- full model: `0d2d032945eeb29d4c19dd0b910c491d8949d53b2bd79074000f1cd3652de716`
- selected fallback: `27fc55fa4e6e39273915df6b2b7b57d043ce3c5abdf1c3eccd9f71db4f2a3842`

A parameter-free historical baseline uses the fraction of earlier-publication
medians already at least 1.0 pChEMBL from the earlier evidence center. It uses no
held-out information and has no fitted coefficients.

## Leakage audit

Every event was checked against the complete source snapshots and event tables
for both earlier phases. The original model ledger contains 27,996 compounds,
41,274 activities, 4,415 assays, 2,736 documents, and five targets. The earlier
independent-validation ledger contains 19,404 compounds, 30,144 activities,
2,802 assays, 1,730 documents, and five targets.

Observed overlap counts were zero for compounds, targets, activities, assays,
and publications against both ledgers. There were also zero duplicate compounds,
zero within-event assay overlaps, and zero within-event publication overlaps.
All historical source measurements precede the held-out publication year. No
cases were excluded as compromised, leaving 94 events from 58 studies and 12
targets.

The full artifact requests 28 molecular or historical-evidence fields. None are
held-out outcomes or post-outcome protocol comparisons. All 12 target categories
are unseen by its one-hot encoder, making this a simultaneous cross-target and
cross-study assessment.

## Prespecified evaluation

The evaluation uses AUROC, PR-AUC, tie-aware Precision@10, Brier score, ECE-10,
and logistic calibration intercept/slope. The tie-aware precision definition
uses the expected number of positives under random ordering at a score tie, so a
constant model does not receive an arbitrary result from CSV row ordering.
Confidence intervals use 2,000 bootstrap replicates of held-out publication
clusters sampled within targets with seed 2026.

### Overall results

| Method | AUROC | PR-AUC | Precision@10 | Brier | ECE-10 | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|---:|---:|---:|
| Frozen full | 0.390 | 0.225 | 0.300 | 0.244 | 0.224 | -1.398 | -0.302 |
| Frozen prevalence | 0.500 | 0.255 | 0.255 | 0.191 | 0.019 | undefined | undefined |
| Historical tenfold fraction | 0.469 | 0.251 | 0.153 | 0.259 | 0.270 | -2.155 | -0.083 |

Frozen-full study-clustered 95% intervals were AUROC 0.240–0.540, PR-AUC
0.136–0.397, Precision@10 0.000–0.500, Brier 0.191–0.307, calibration intercept
-2.633 to -0.568, and calibration slope -1.029 to 0.221.

Compared with the frozen prevalence fallback, the full model changed AUROC by
-0.110 (95% CI -0.265 to 0.043), PR-AUC by -0.030 (-0.070 to 0.082), and Brier
by +0.054 (+0.015 to +0.092; positive is worse). Its apparent Precision@10 gain
of +0.045 had a wide interval (-0.267 to +0.260). There is no defensible overall
improvement.

### Natural products

The subgroup includes 20 compounds from 17 held-out studies, with seven
disagreements. Results remain underpowered and are not representative of
natural products generally.

| Method | AUROC | PR-AUC | Precision@10 | Brier |
|---|---:|---:|---:|---:|
| Frozen full | 0.209 | 0.261 | 0.100 | 0.346 |
| Frozen prevalence | 0.500 | 0.350 | 0.350 | 0.241 |
| Historical tenfold fraction | 0.407 | 0.336 | 0.300 | 0.357 |

Frozen-full intervals were AUROC 0.107–0.338, PR-AUC 0.175–0.383,
Precision@10 0.000–0.300, and Brier 0.243–0.444. Calibration slope was negative
(-1.351; 95% CI -2.520 to -0.627), indicating inverse ordering in this small
subgroup rather than transferable discrimination.

At the fixed, non-tuned probability threshold of 0.5, all seven natural-product
disagreements were missed and one negative was flagged. Only one disagreement,
`CHEMBL426110` against `CHEMBL1951`, entered the full model's top ten, at rank 8
with predicted risk 0.258. The other six disagreement cases ranked 12–20.

## Conclusion

This evaluation is a clear negative result. The frozen full model does not
generalize reliably to the expanded unseen-target cohort and is especially poor
for the natural-product subgroup. The constant fallback is not discriminative;
its better Brier score reflects conservative calibration around the development
prevalence, not useful candidate ranking. No model is selected as a winner and
no parameter or threshold is changed.

The observed 25.5% overall and 35.0% natural-product disagreement fractions are
properties of an enriched repeated-measurement sample. They must not be treated
as real-world assay-disagreement prevalence. Disagreement itself is an assay
reproducibility outcome, not evidence of efficacy, safety, biological failure,
or clinical utility.

## Reproduction

From `backend`:

```powershell
python -m app.expanded_validation --events ..\experiments\assay_disagreement_expansion\expansion_events.csv --measurements ..\experiments\assay_disagreement_expansion\used_measurements.csv --protocol ..\experiments\expanded_validation_protocol.json --development-specs ..\experiments\multitarget_targets.json --development-events ..\experiments\assay_disagreement_generalization\dataset.csv --prior-validation-specs ..\experiments\independent_validation_targets.json --prior-validation-events ..\experiments\independent_assay_validation\independent_dataset.csv --full-model ..\experiments\assay_disagreement\assay_disagreement_model.joblib --selected-model ..\experiments\assay_disagreement_generalization\study_generalized_model.joblib --output-dir ..\experiments\expanded_assay_validation
python -m pytest -q --basetemp .pytest-tmp-expanded-validation
```

The output manifest records all input and output checksums. Prediction rows,
metrics, confidence intervals, leakage exclusions, natural-product failures,
comparison plots, and the generated report are preserved alongside it.
