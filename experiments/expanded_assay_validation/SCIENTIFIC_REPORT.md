# Expanded independent validation of frozen assay-disagreement models

## Independence audit

All **94 eligible events** passed compound, target, activity, assay, publication, temporal, and within-event leakage checks against both the original model corpus and the earlier independent validation. The audit found zero compromised events. Model SHA-256 values were identical before and after inference.

The frozen full artifact uses only molecular descriptors and historical evidence fields. The selected study-generalized artifact remains the constant development prevalence. No fitting, recalibration, threshold tuning, or winner selection was performed.

## Dataset

The locked expansion contains 94 compounds, 24 disagreements, 58 held-out studies, and 12 unseen targets. The observed disagreement fraction (0.255) is a property of this enriched repeated-measurement sample and is **not** an estimate of real-world prevalence.

## Overall performance

| Method | AUROC | PR-AUC | Precision@10 | Brier | ECE-10 | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|---:|---:|---:|
| Frozen full | 0.390 | 0.225 | 0.300 | 0.244 | 0.224 | -1.398 | -0.302 |
| Frozen prevalence | 0.500 | 0.255 | 0.255 | 0.191 | 0.019 | undefined | undefined |
| Historical tenfold fraction | 0.469 | 0.251 | 0.153 | 0.259 | 0.270 | -2.155 | -0.083 |

Frozen-full study-clustered 95% intervals: `{"auroc": [0.2403670634920635, 0.5404080478916179], "brier": [0.19146074616603484, 0.30749957204541406], "calibration_intercept": [-2.632708375281628, -0.5677646702667147], "calibration_slope": [-1.0288127122537232, 0.2208486539690419], "ece_10": [0.17411165144668939, 0.33139327872224233], "log_loss": [0.5800726556739062, 0.8695723567193451], "positives": [16.0, 33.0], "pr_auc": [0.13577261425949685, 0.39672588544309284], "precision_at_10": [0.0, 0.5], "prevalence": [0.16666666666666666, 0.36363636363636365]}`.

## Natural products

The subgroup has 20 compounds and 7 disagreements. Frozen-full AUROC is `0.20879120879120877`, PR-AUC `0.26051305137439595`, Precision@10 `0.1`, and Brier `0.34582609279389753`. Study-clustered intervals are `{"auroc": [0.10714285714285715, 0.3384615384615385], "brier": [0.24265098781707742, 0.44440706214717884], "calibration_intercept": [-5.075604144086427, -1.6779814388713141], "calibration_slope": [-2.51978381928344, -0.627377815915854], "ece_10": [0.2707088902628231, 0.49279533523614155], "log_loss": [0.717726226891265, 1.3395146800714681], "positives": [4.0, 11.0], "pr_auc": [0.17536630036630035, 0.3829153154537769], "precision_at_10": [0.0, 0.3], "prevalence": [0.2222222222222222, 0.4583333333333333]}`.

At the fixed, non-tuned 0.5 descriptive threshold, classifications are `{"false_positive": 1, "missed_disagreement": 7, "true_negative": 12}`. The accompanying error table preserves every natural-product compound, score, rank, outcome, target, publication, and potency shift.

## Conclusion

The expanded evaluation does not establish reliable generalization: uncertainty includes weak or non-discriminatory performance under the prespecified criteria. The negative or inconclusive result is retained without selecting or modifying a model.

No method is selected as a winner from this evaluation. Disagreement is an assay-reproducibility outcome, not evidence of efficacy, safety, biological failure, or clinical utility.

## Limitations

- The test set was deliberately enriched for independently repeated measurements and is not prevalence-representative.
- All 12 targets are unseen by the frozen target encoder, so this is simultaneous cross-target and cross-study validation.
- Twenty natural products and seven subgroup positives yield wide, unstable confidence intervals.
- Publication identifiers proxy independent studies; laboratory and full protocol metadata remain incomplete.
- Calibration estimates are descriptive and no recalibration was permitted.
