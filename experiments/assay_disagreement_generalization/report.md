# Study-grouped assay-disagreement generalization analysis

## Decision protocol

The 2020+ evaluation was already summarized in the preceding experiment, so it is locked but not analyst-naive. No feature, candidate, hyperparameter, calibration mapping, or selection rule in this phase uses its outcomes. Selection uses only leave-one-study-out predictions from publications through 2019, minimizing equal-study-weighted Brier score. Natural-product status is never a predictor.

The prespecified outcome remains a ≥1.0-log10 difference between earlier and later exact-IC50 medians from document- and assay-disjoint records in one compatible context. It is an assay disagreement signal, not drug efficacy, safety, biological failure, or clinical outcome.

## Development selection

There were 212 development compounds from 141 held-out studies; 112 studies contained one compound. The selected method was `constant_prevalence` by minimum equal-study-weighted Brier score from leave-one-study-out development predictions. All candidate LOSO metrics are in `development_model_metrics.csv`.

## Locked evaluation

| Model | AUROC | 95% study-cluster CI | PR-AUC | 95% study-cluster CI | Brier | Calibration slope |
|---|---:|---:|---:|---:|---:|---:|
| Current full logistic | 0.642 | [0.414858059445303, 0.6820809653248475] | 0.609 | [0.3114092209229677, 0.6528511544373704] | 0.296 | 0.258 |
| Development-selected `constant_prevalence` | 0.500 | [0.5, 0.5] | 0.530 | [0.29223270440251575, 0.6594710892768786] | 0.336 | not identifiable |

The paired selected-minus-current differences and clustered intervals are recorded in `results.json`. The selected model's status is **not_supported_for_reliable_prediction**: At least one prespecified study-generalization or calibration condition failed.

## Study bias

Binary disagreement ICC was 0.2883610901772933. CHEMBL5729139 contained 131 compounds and 102 disagreements (observed rate 0.779); mean selected risk was 0.236, while the current full model's mean risk was 0.549. The current full model's AUROC/PR-AUC after excluding this study were 0.4826692981374121/0.3580382113957865. This study is a protocol/batch shift the model cannot identify from earlier compound evidence alone.

## Natural products

The evaluation subgroup contained 60 structure-verified natural products and 19 disagreements. The prior full model had AUROC 0.48267008985879334 and PR-AUC 0.33626511367971135; the development-selected fallback had AUROC 0.5 (95% study-cluster CI [0.5, 0.5]), PR-AUC 0.31666666666666665 (CI [0.2096078976640712, 0.4406779661016949]), and Brier score 0.22292037498516673. 4 compounds were outside the development descriptor domain. The dominant overall EGFR study contained no verified natural products, while NP compounds had materially longer and more heterogeneous evidence histories. Compound errors and descriptor shifts are in dedicated CSVs.

## Interpretation

Study-grouped selection did not establish reliable generalization. Study identity and protocol shift remain stronger determinants of disagreement than transferable compound-history signals.

## Limitations

- Only 212 development compounds and 50 development disagreements were available after strict duplicate filtering; 112/141 development studies were singletons.
- Publication ID is an imperfect proxy for laboratory/protocol independence, and structured substrate, construct, detection, incubation, and laboratory variables remain sparse.
- The locked evaluation was observed in the previous phase, so this is not a fully analyst-blinded replication.
- One EGFR publication dominates positive evaluation outcomes; molecule-level confidence intervals would be anticonservative.
- Natural-product estimates remain small and target-confounded. NPASS/COCONUT annotations are subgroup labels only.
- Cross-organism generalization is not estimable from these five human-target contexts.
