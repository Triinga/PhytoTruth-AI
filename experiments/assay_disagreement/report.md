# Evidence Transferability and Assay Disagreement Prediction

## Question and outcome

Can information already available for a compound predict whether an **independent, later ChEMBL study** will differ materially from its earlier exact-IC50 evidence? The primary disagreement outcome is an absolute shift of at least 1.0 log10 pActivity (ten-fold potency). This is an assay-reproducibility signal, not evidence of biological failure, efficacy, safety, or clinical outcome.

Only exact, normalized IC50 records in each target's prespecified compatible human single-protein context were used. Earlier and held-out activity, document, and assay identifiers are disjoint. Units, endpoint, target organism, assay type, and BAO format are therefore not pooled across incompatible contexts.

## Design

- Development outcomes: latest independent publication through 2019; n=212, disagreements=50.
- Untouched evaluation outcomes: latest independent publication from 2020; n=319, disagreements=169.
- Molecules are disjoint across development and evaluation. Features use only structure and evidence preceding that molecule's held-out publication. Held-out assay metadata is never a predictor.
- The model is a prespecified L2 logistic regression. Comparators and ablations were fixed before evaluation. No evaluation-set tuning was performed.
- Intervals resample held-out study clusters within target; molecules from the same publication are not treated as independent studies.

| Target | Development | Evaluation | Evaluation disagreements | Evaluation NPs |
|---|---:|---:|---:|---:|
| CHEMBL220 | 49 | 42 | 11 | 20 |
| CHEMBL203 | 81 | 212 | 132 | 12 |
| CHEMBL230 | 42 | 25 | 8 | 16 |
| CHEMBL284 | 16 | 25 | 16 | 7 |
| CHEMBL335 | 24 | 15 | 2 | 5 |

## Independent evaluation

The constant development-prevalence baseline obtained AUROC 0.5 and PR-AUC 0.530. The full model obtained AUROC 0.642 (study-clustered 95% CI [0.41350239920099013, 0.6864486163132736]) and PR-AUC 0.609 (95% CI [0.3032084700239042, 0.6539707547838226]). Its Brier score was 0.296, calibration intercept 0.455, and calibration slope 0.258.

The verified-natural-product subgroup had n=60 and 19 disagreements; full-model AUROC=0.48267008985879334 and PR-AUC=0.33626511367971135. This subgroup is exploratory because study and positive-event counts are small.

**Validation decision:** `not_supported_for_reliable_prediction`. The prespecified support rule required the study-clustered AUROC interval to exclude 0.5 and the PR-AUC interval to exceed evaluation prevalence. At least one criterion failed, so the serialized model is an auditable research artifact and must not be presented as a validated predictor.

Full numerical baselines, ablations, target-specific results, calibration, and uncertainty are in `model_metrics.csv` and `results.json`. `predictions.csv` preserves event provenance and per-compound bootstrap risk intervals. Coefficients are standardized associations, not causal effects.

## Transferability findings

The prespecified full model did not establish reliable prediction beyond simple baselines. Full PR-AUC=0.609; target-only=0.515; molecular-only=0.512; earlier-evidence-only=0.559. Performance is reported as predictive association only; no evidence component is claimed transferable unless its ablation improves independently evaluated performance with compatible uncertainty.

The main source of uncertainty is study-level protocol shift: CHEMBL5729139 is the largest held-out publication cluster (131 compounds, 102 disagreements; 60.4% of all evaluation disagreements). A model that cannot see a future protocol cannot fully anticipate such a batch shift. The sensitivity outcome crossing pActivity 6.0 is reported descriptively and is not substituted for the prespecified ten-fold outcome.

## Limitations

- ChEMBL assay descriptions do not provide uniformly structured protocol variables such as substrate concentration, incubation time, construct, detection technology, or laboratory; residual protocol heterogeneity remains.
- The split is publication-temporal but uses the current ChEMBL and NPASS/COCONUT releases. Natural-product status is used only for subgroup evaluation, never as a predictor.
- Latest-publication events are retrospective. Database deposition and duplicated literature values can still create selection effects despite document/assay disjointness and duplicate flags.
- Only five human targets and one compatible exact-IC50 context per target were studied. Cross-organism transferability is not estimable.
- Natural-product estimates are underpowered and cannot establish separate performance.
- Disagreement does not identify which study is correct and must not be interpreted as clinical failure.
