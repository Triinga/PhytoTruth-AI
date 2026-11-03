# Independent assay-disagreement validation

## Independence contract

This analyst-independent target replication used **5 previously unused ChEMBL targets** and held-out publications from 2020 onward. Every compound, activity, assay, and publication appearing anywhere in the previous five normalized target snapshots was excluded. The frozen model artifacts were loaded without fitting, calibration, threshold selection, or parameter changes; their SHA-256 values were verified unchanged after scoring.

All validation targets are outside the frozen model's five-target training domain. The existing one-hot encoder therefore treats target identity as unknown; results measure cross-target as well as cross-study generalization.

## Dataset

The final dataset contains **90 compounds**, **7 disagreements**, and **37 held-out publications**. It includes **5 structure-verified natural products** with **0 disagreements across 5 studies**. These events are backed by **277 traceable measurements** (180 prior and 97 held out).

| Target | Name | Events | Disagreements | Studies | Natural products | NP disagreements |
|---|---|---:|---:|---:|---:|---:|
| CHEMBL1824 | Receptor tyrosine-protein kinase erbB-2 | 2 | 0 | 2 | 0 | 0 |
| CHEMBL1862 | Tyrosine-protein kinase ABL1 | 8 | 0 | 3 | 0 | 0 |
| CHEMBL1871 | Androgen receptor | 16 | 3 | 12 | 1 | 0 |
| CHEMBL206 | Estrogen receptor | 9 | 0 | 8 | 3 | 0 |
| CHEMBL279 | Vascular endothelial growth factor receptor 2 | 55 | 4 | 12 | 1 | 0 |

## Frozen-model performance

| Model | AUROC | Study-clustered 95% CI | PR-AUC | Study-clustered 95% CI | Brier | ECE-10 | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Frozen full logistic | 0.4939759036144578 | [0.37034908026755853, 0.8555875576036868] | 0.09345500232477297 | [0.0771926406926407, 0.3333333333333333] | 0.42075286497359043 | 0.5377740495826734 | -2.447976775976694 | -0.06841221286352502 |
| Frozen development-selected fallback | 0.5 | [0.5, 0.5] | 0.07777777777777778 | [0.03333333333333333, 0.13759698275862065] | 0.0967149242514141 | 0.1580712788259958 | None | None |

## Natural-product subgroup

The frozen full model had AUROC `None`, PR-AUC `None`, Brier `0.09653298855938017`, and study-clustered intervals `{'brier': [0.020169321576549878, 0.14092841855788113], 'log_loss': [0.1115099813300775, 0.45764269260157076], 'ece_10': [0.0995387760159841, 0.3575551908135688]}`. The prespecified subgroup sufficiency criteria are reported in `results.json`; metrics remain descriptive if those criteria are not met.

## Decision

**not_independently_supported**. Criteria: `{"auroc_ci_lower_above_half": false, "brier_no_worse_than_frozen_selected": false, "minimum_dataset_met": false, "pr_auc_ci_lower_above_prevalence": false}`.

The prespecified overall minimum was not met by `{'events': 10, 'disagreements': 13, 'heldout_studies': 0}`. Natural-product deficits were `{'events': 25, 'disagreements': 10, 'heldout_studies': 0}`. These deficits describe the minimum additional eligible observations needed for the planned analysis; they do not guarantee adequate statistical power.

This result concerns prediction of later-assay disagreement only. It does not identify which assay is correct and does not imply drug efficacy, safety, biological failure, or clinical performance.

## Limitations

- The targets were frozen without inspecting disagreement outcomes, but chosen partly for ChEMBL activity volume and biological relevance; they are not a random sample of pharmacological targets.
- All targets are unseen by the frozen target encoder, so cross-target shift is deliberately severe and target-specific intercepts cannot transfer.
- Publication identifiers are imperfect proxies for independent laboratories and protocols.
- Natural-product membership is structure-verified NPASS/COCONUT membership, not proof of natural origin for the specific assayed sample.
- The latest-publication construction is retrospective and ChEMBL deposition can be incomplete or delayed.
- Confidence intervals preserve held-out-study clustering within target; few clusters or positives can still make them unstable.

## Reproduction

```powershell
python -m app.pipeline independent-assay-validation --specs ../experiments/independent_validation_targets.json --protocol ../experiments/independent_validation_protocol.json --previous-specs ../experiments/multitarget_targets.json --previous-events ../experiments/assay_disagreement_generalization/dataset.csv --full-model ../experiments/assay_disagreement/assay_disagreement_model.joblib --selected-model ../experiments/assay_disagreement_generalization/study_generalized_model.joblib --output-dir ../experiments/independent_assay_validation --bootstrap-iterations 2000 --seed 2026
```
