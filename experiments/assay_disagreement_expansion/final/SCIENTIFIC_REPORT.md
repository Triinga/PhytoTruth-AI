# Assay Disagreement Dataset Expansion

Generated 2026-10-21 under `assay-disagreement-dataset-expansion-v1`.

## Cohort

The strictly independent expansion contains **94 compound-target events** across **58 held-out publications** and **12 targets**. There are **24 disagreements** (absolute prior-versus-independent median shift >= 1.0 pChEMBL) and **70 negatives**. Exact standardized-structure matches identify **20 verified natural products**, including **7 disagreements**.

## Scientific interpretation

All labels use exact IC50 measurements within one target and a frozen compatible assay context. Historical and held-out evidence have disjoint ChEMBL assay and document identifiers. Incompatible endpoints, censored relations, normalization mismatches, potential duplicates, and records carrying ChEMBL validity comments are excluded rather than treated as disagreements.

Protocol comparisons are post-outcome explanatory annotations. A documented difference is a possible explanation only; it is not causal evidence. `unexplained_with_partial_documented_similarity` means at least four reported dimensions agree while other protocol fields may remain unknown. `unresolved_due_to_incomplete_protocol_metadata` means ChEMBL does not report enough protocol detail to interpret the shift. Missing fields are never inferred.

## Positive-case protocol classification

- `documented_protocol_difference_possible_explanation`: 4
- `unexplained_with_partial_documented_similarity`: 3
- `unresolved_due_to_incomplete_protocol_metadata`: 17

## Label quarantine

These new labels are reserved as an independent dataset. They were not used to train, tune, select, or evaluate a model in this phase. No model artifact is loaded by this pipeline.

## Limitations

- ChEMBL assay descriptions are often summaries, not full experimental protocols.
- Structural membership in NPASS/COCONUT verifies database membership, not biosynthetic origin or clinical efficacy.
- Latest-year evidence can include several records from one publication and is therefore one study-level event, not several independent confirmations.
- A one-log-unit threshold is a prespecified operational definition of disagreement, not proof of biological failure.
