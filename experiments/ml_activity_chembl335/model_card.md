# PhytoTruth Molecular Bioactivity Model

Generated: 2026-10-20T12:19:12.510239+00:00

## Outcome

The selected learned model is **xgboost**. It outperformed the mean baseline on scaffold-held-out MAE.

This model predicts experimental pActivity for the defined ChEMBL assay context. It does not predict clinical efficacy, safety, or a probability that a compound is a drug.

## Dataset

- Target/context: `bee3486ec6fbd9a6ba15` — {"assay_type": "B", "bao_format": "BAO_0000357", "standard_type": "IC50", "target_chembl_id": "CHEMBL335", "target_organism": "Homo sapiens"}.
- 2915 exact compatible measurements aggregated to 2444 standardized structures.
- 362 structures are verified natural products through the supplied structure crossmatch (NPASS for this executed run).
- 1086 Bemis–Murcko scaffold groups; 5 outer folds.
- Replicate measurements were aggregated by standard InChIKey before splitting, preventing the same standardized structure from crossing folds.

## Scaffold-held-out results

| Model | MAE | RMSE | R² | Natural-product MAE | 90% interval coverage | Outside AD |
|---|---:|---:|---:|---:|---:|---:|
| mean_baseline | 0.738 | 0.925 | -0.006 | 0.662 | 0.900 | 279 |
| random_forest | 0.455 | 0.597 | 0.581 | 0.396 | 0.905 | 306 |
| xgboost | 0.432 | 0.575 | 0.612 | 0.368 | 0.924 | 315 |

Selected-model MAE 95% scaffold-cluster bootstrap CI: [0.408, 0.461].

## Validation and uncertainty

- Outer test folds contain scaffolds absent from their model-fitting and calibration partitions.
- Each outer training set is split again by scaffold into fitting and uncertainty-calibration partitions.
- The 90% interval is based on calibration absolute residuals and inflated when Tanimoto similarity falls below the fold-specific applicability threshold.
- Applicability thresholds are the fifth percentile of nearest-neighbor Morgan similarity within training data. This is an operational heuristic, not proof that an in-domain prediction is correct.
- Selected-model empirical 90% interval coverage was 0.924; mean interval width was 2.069 pActivity units.

## Model selection

The learned model with the lowest outer-fold MAE was selected. Mean-baseline MAE was 0.738; selected-model MAE was 0.432. No hyperparameter search used outer-fold outcomes.

XGBoost-minus-Random-Forest MAE was -0.023 (95% paired scaffold-cluster bootstrap CI [-0.038, -0.010]).

The interval-width versus absolute-error Spearman diagnostic was -0.003. Near-zero association means the per-molecule width did not meaningfully rank errors, even though aggregate coverage was close to its target.

## Limitations

- Model selection and performance estimation use the same fixed cross-validation run; an untouched external target-specific test set is not available.
- ChEMBL measurements are heterogeneous literature evidence despite strict metadata filtering and per-structure aggregation.
- Scaffold splits are deliberately difficult but do not reproduce every real deployment distribution shift.
- Prediction intervals are cross-validated residual diagnostics, not guaranteed prospective coverage, especially outside the applicability domain.
- Natural-product subgroup metrics are observational; natural-product status is not used as a predictive feature.
- The Evidence Reversal Experiment remains inconclusive/negative and is not evidence that its reliability heuristic is validated.

## Reproduction

Use `command.txt`, `manifest.json`, `environment.json`, and `feature_config.json` in this directory. The serialized model is `best_model.joblib`.
