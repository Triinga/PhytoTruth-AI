# Molecular Bioactivity Model Experiment

## Scope

This experiment predicts pActivity for one compatible ChEMBL assay context. It does not predict clinical efficacy, safety, ethnobotanical validity, or probability of drug success. The Evidence Reversal reliability heuristic is not used as a validated feature because its experiment was inconclusive/negative.

## Dataset construction

1. Read normalized real ChEMBL JSONL.
2. Keep exact, finite measurements from a single context: target `CHEMBL335`, human IC50, binding assay, BAO single-protein format.
3. Standardize each structure using the pipeline's RDKit parent/uncharging procedure.
4. Aggregate replicate activity to the median for each standard InChIKey before splitting.
5. Attach NPASS/COCONUT membership from the structure-verified crossmatch. Membership is used only for subgroup evaluation, never as a predictive feature.

The executed dataset contained 2,915 measurements, 2,444 unique standardized structures, 1,086 Bemis–Murcko scaffolds, and 362 NPASS-verified natural products. Nineteen selected records lacked a usable structure and were excluded. Other assay contexts were not combined.

## Features

- Morgan fingerprint, radius 2, 2,048 bits.
- RDKit molecular weight, logP, TPSA, hydrogen-bond donors and acceptors, rotatable bonds, fraction Csp3, ring count, heavy-atom count, and formal charge.
- Natural-product annotations, evidence reliability scores, activity years, compound IDs, and source identifiers are not model inputs.

## Validation

Five-fold Bemis–Murcko scaffold validation keeps every standardized structure and scaffold in exactly one outer test fold. Inside each outer training partition, a second scaffold-grouped 80/20 split separates model fitting from uncertainty calibration. The same splits and fixed hyperparameters are used for the mean baseline, Random Forest, and XGBoost.

No hyperparameter search uses outer-fold outcomes. XGBoost is selected only because it has the lowest outer-fold MAE among the learned models. Its paired advantage over Random Forest is reported with a scaffold-cluster bootstrap rather than inferred from point estimates alone.

## Uncertainty and applicability domain

Each outer fold calculates a finite-sample 90% absolute-residual conformal quantile on its calibration scaffolds. Test intervals use that half-width and are widened below the applicability threshold.

The applicability threshold is the fifth percentile of nearest-neighbor Morgan Tanimoto similarity among fitting structures. A prediction below the threshold is marked outside domain. This flag is an operational distance diagnostic—not proof that an in-domain prediction is reliable. The serialized final model stores packed training fingerprints, the threshold, and the cross-validated residual half-width.

## Executed results

| Model | MAE | RMSE | R² | NPASS subgroup MAE | 90% coverage |
|---|---:|---:|---:|---:|---:|
| Mean baseline | 0.738 | 0.925 | -0.006 | 0.662 | 0.900 |
| Random Forest | 0.455 | 0.597 | 0.581 | 0.396 | 0.905 |
| XGBoost | 0.432 | 0.575 | 0.612 | 0.368 | 0.924 |

XGBoost was best on the prespecified primary metric. Its MAE 95% scaffold-cluster interval was 0.408–0.461. XGBoost-minus-Random-Forest MAE was -0.023 with a paired scaffold-cluster interval of -0.038 to -0.010.

For XGBoost, in-domain MAE was 0.420 and outside-domain MAE was 0.511; 315 of 2,444 outer-fold predictions were outside domain. The mean 90% interval width was 2.069 pActivity units. Interval width had essentially no rank association with absolute error (Spearman -0.003), so aggregate coverage was acceptable but molecule-specific uncertainty ranking was not informative. This negative calibration diagnostic must remain visible.

## Reproduction

```powershell
cd backend
python -m app.pipeline train-activity `
  --activities ../data/processed/chembl335_ic50.normalized.jsonl `
  --natural-products ../data/processed/chembl335_ic50.natural-products.jsonl `
  --output-dir ../experiments/ml_activity_chembl335 `
  --folds 5 --bootstrap-iterations 1000 --seed 2026
```

The output contains the serialized model, dataset/fold assignments, out-of-fold predictions, fold metrics, feature configuration, environment versions, input/model checksums, comparison SVG, and model card.

## Limitations

- There is no independent external target-specific test set.
- Selecting between two learned models on the same cross-validation estimates can introduce selection optimism.
- Assay metadata filtering reduces but cannot eliminate protocol heterogeneity.
- Natural-product subgroup performance is observational and target-specific.
- Split-conformal diagnostics are retrospective; future chemical series may shift beyond the represented domain.
- Predictions outside the domain must not be used without expert review and new experiments.
