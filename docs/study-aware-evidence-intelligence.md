# Study-aware evidence intelligence

## Purpose

This module turns a point prediction into an auditable research-planning decision. It estimates a future-study pActivity distribution, the probability of a ten-fold or larger shift from historical evidence, when the system should abstain, and which independent experiment is expected to reduce uncertainty most. None of these outputs are efficacy, safety, or clinical claims.

## Source and selection

The compiler consumes the executed PhytoTRUST-M1 development-only artifact. It compares six prespecified methods using two existing validation schemes: leave-one-publication-out study-grouped validation and strictly forward-temporal validation. Selection minimizes the worst negative-log-likelihood rank across both schemes, then mean rank and calibration error. This selects Student-t regression: rank 1 for study-grouped NLL and rank 2 for forward-temporal NLL.

This is a secondary retrospective analysis of an already inspected development dataset. It is not independent confirmation.

## Disagreement probability

For every forward-temporal prediction, the fitted Student-t predictive distribution is integrated outside historical median pActivity ±1.0. This yields a probability of a ten-fold or larger future-study shift. The outcome label is used only for retrospective evaluation. It is never read by recommendation prioritization.

AUROC, average precision, and Brier score are accompanied by publication-cluster bootstrap intervals. The executed sample contains only 25 events in eight publications and three disagreements. Consequently, the AUROC point estimate of 0.826 has a wide 95% study-cluster interval of 0.375–0.958 and is not evidence of a validated classifier.

## Abstention

Predictions are ordered by 90% interval width, using no held-out outcome. Retrospective MAE and coverage are then reported at 100%, 80%, 60%, and 40% retention. At 60% retention, MAE decreases from 0.517 to 0.317 and observed 90% coverage rises from 84% to 100%. This promising selective-prediction pattern requires new independent replication.

## Expected information gain

The planning score approximates entropy reduction from one independent observation under a Gaussian latent-mean approximation to the Student-t interval. It increases with predictive uncertainty and decreases with the number of prior independent studies. The recommendation uses only prediction-time fields: the prediction, interval, historical median, prior-study count, and historical protocol missingness.

This utility prioritizes information collection. It is not a probability that an experiment will succeed, and it has not been prospectively validated.

## Reproduction

Run `experiments/evidence_intelligence/command.txt` from `backend/`. The output includes evaluated predictions, recommendations, results, source checksums, the explicit outcome-blind field contract, and the deployment guardrail.
