# Architecture

## System boundary

PhytoTruth is a research-prioritization system. It does not infer clinical efficacy, diagnose disease, or replace confirmatory experiments.

```text
COCONUT CSV/SDF ─┐
NPASS exports ───┼─> immutable raw zone ─> canonical structures ─┐
ChEMBL REST API ─┘                                               │
                                                                v
                                                     evidence normalization
                                                  (units, relations, provenance)
                                                                │
                            ┌───────────────────────────────────┴────────────┐
                            v                                                v
                  RDKit Morgan fingerprints                       reliability engine
                  + Murcko scaffold groups                    dispersion / conflicts /
                            │                                  organisms / duplicates
                            v                                                │
               Random Forest / XGBoost                                  v
                  scaffold-grouped CV                            evidence-aware ranking
                            └───────────────────┬────────────────────────────┘
                                                v
                                      FastAPI provenance API
                                                v
                                      React research workspace
```

## Data contracts

Raw records are append-only and retain source, release/retrieval date, source identifier, assay identifier, document identifier, original units, relation and validity flags. Processed records add canonical SMILES and normalized nanomolar values without deleting original fields.

An experimental result has `kind=measured`; a model output has `kind=predicted`, model version, training snapshot and applicability-domain metadata. The two are never rendered with the same visual label.

## Reliability model

The MVP reliability score is deterministic and inspectable. It increases with independent assay and document support; it decreases with median absolute deviation, threshold-crossing contradictions, possible duplicates, source validity comments and mixed organisms. The evidence ranking uses median potency, not the best reported value, and includes explicit uncertainty and support terms.

This heuristic is intentionally not presented as a calibrated probability. A later calibration stage can learn weights against prospective replication outcomes.

## Leakage controls

- Molecular records are grouped before splitting, so replicate measurements cannot appear on both sides of a model-validation fold.
- ML validation aggregates standardized structures before using Bemis–Murcko scaffold groups; fitting, uncertainty calibration and outer testing use disjoint scaffold sets.
- Evidence-reversal evaluation freezes rankings from training evidence and scores only held-out consensus.
- When dates are available, the preferred experiment is a temporal split by publication date.
- Hyperparameter tuning, if added, must be nested inside the training partition.

## Evaluation

Regression: MAE, RMSE and R² from out-of-fold predictions, with scaffold-cluster bootstrap confidence intervals. Nested scaffold calibration estimates residual intervals. Morgan-Tanimoto nearest-neighbor distance marks predictions outside the applicability domain; the flag is a distance heuristic, not calibrated confidence.

Ranking reversal: nDCG@k, enrichment factor at k, precision@k and Kendall τ against held-out median activity. Report molecule count, active threshold, split policy and confidence intervals with every result.

## Extensibility

The ingestion adapters stream bulk files, allowing large COCONUT exports without holding them in memory. New sources implement the same provenance contract. A production deployment should add object storage for immutable releases, Postgres for normalized evidence, an experiment registry (for example MLflow), and background workers for fingerprinting.

