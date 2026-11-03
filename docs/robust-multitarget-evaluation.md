# Robust multi-target ranking evaluation

## Question and claim boundary

This benchmark asks whether reliability-aware historical evidence improves ranking of molecules against later published IC50 evidence. It does not test efficacy, safety, clinical utility, or biological truth. A lower p-value is never an optimization target.

## Frozen data design

The runner discovers local target files only when both normalized ChEMBL measurements and natural-product crossmatches exist. Eligibility depends only on coverage: at least 25 molecules must have compatible measurements in multiple publication years. Whole targets are assigned deterministically to development (60%) or evaluation (40%) from the target ID and seed, without examining outcomes.

For each molecule, all measurements in the latest publication year are held out. Every model input is calculated from strictly earlier evidence. Activity, assay, and document identifiers must be disjoint across the historical and held-out sets. Natural-product membership is descriptive and is not a model feature.

## Development-only optimization

The candidate predicts a residual adjustment to historical median potency with ridge regression. Candidate feature families and ridge strengths are compared by leave-one-development-target-out validation. Molecules are weighted so each target has equal total training weight. The frozen development gate requires both nDCG@10 and Kendall tau to have a target-bootstrap confidence interval above zero and Holm-adjusted exact target sign-flip p<0.05. Failure retains median potency as the deployment guardrail.

## Locked evaluation and inference

The selected configuration is fit on development targets, then evaluated once on untouched targets. Metrics are calculated within target and macro-averaged. The independent unit for inference is the biological target, not the molecule:

- 95% intervals resample target clusters.
- Two-sided p-values use exact sign flips of target-level paired effects.
- Holm correction covers the four reported ranking metrics.
- Paired Hedges g and target wins/ties/losses accompany p-values.

With eight evaluation targets, the minimum attainable two-sided exact p-value is 2/2^8 = 0.0078125. Consequently, this run cannot establish p<=0.001. At least 11 nonzero independent target effects are needed even for that p-value to be mathematically attainable; substantially more may be needed for adequate power at realistic effect sizes.

## Ablations

Support, consistency, reliability/uncertainty, and recent-drift feature groups are removed from the frozen selected model one group at a time. Each full-minus-ablation effect uses the same target-level inference. A group absent from the selected feature family is marked not applicable instead of being reported as a zero-effect experiment.

## Reproduction

Run `experiments/robust_multitarget_evaluation/command.txt` from `backend/`. The output directory contains input SHA-256 hashes, the deterministic target partition, model coefficients, target-level effects, ablation inference, JSON results, and a generated report. This is a retrospective benchmark over reused public database targets, not a prospective or independent-laboratory confirmation.
