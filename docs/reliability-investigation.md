# Evidence Reliability Engine validation protocol

This investigation asks whether the existing reliability adjustments improve prediction of later, document-disjoint ChEMBL bioactivity evidence beyond robust median potency. It does not evaluate drug efficacy or clinical safety.

## Frozen design

1. Select one compatible target, organism, activity type, assay type, and BAO format without inspecting activity values.
2. For every eligible molecule, reserve all measurements from its latest publication year. Earlier publications form the scoring evidence; publication documents cannot appear on both sides.
3. Stratify molecules only by NPASS verification and deterministically hash molecule IDs into 60% development and 40% evaluation sets.
4. On development data only, compare the existing PhytoTruth score with median potency. The evidence adjustment passes only if nDCG@10 and Kendall tau-b both have positive paired 95% interval lower bounds and two-sided score-swap p-values below 0.05.
5. Freeze the choice, then report the untouched-for-revision evaluation results. Because the parent cohort was used in an earlier experiment, this is an internal revision holdout—not an external validation dataset.

## Component tests

Each component is assessed against two future-evidence outcomes:

- Signed residual (`held-out median - training median`) tests whether a score adjustment points in the correct direction.
- Absolute residual tests whether a risk signal identifies future prediction error.

Spearman correlations receive molecule-bootstrap 95% intervals and permutation p-values. These associations do not convert heuristic scores into probabilities.

## Error attribution

A harmful reversal is a molecule pair ordered correctly by training median potency but incorrectly after the current PhytoTruth adjustment. Attribution uses the exact additive score terms: assay support, explicit inconsistency, composite reliability, and bootstrap uncertainty. The reported dominant component is the largest pairwise term opposing the original median ordering; it is a deterministic decomposition, not a causal claim.

## Interpretation boundary

If the development gate fails, no new weights are fitted. Median potency becomes the defensible ranking guardrail while reliability fields remain descriptive warnings. Establishing a learned adjustment requires substantially more natural products, multiple targets, richer protocol and laboratory metadata, and an external publication/laboratory holdout.
