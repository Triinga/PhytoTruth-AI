# Ethnopharmacological convergence: traditional antidiabetic use vs predicted PTP1B activity

**Protocol:** `ethno-convergence-v1` — frozen before any outcome was computed.

## Verdict: NOT SUPPORTED

> Do plants with a documented traditional antidiabetic use yield natural products with higher model-predicted PTP1B activity than Duke-covered control plants without such a use?

Plants with a documented traditional antidiabetic use had a median predicted pChEMBL **0.003 log units lower** than Duke-covered control plants.

| Quantity | Value |
|---|---:|
| Traditional antidiabetic plants | 121 |
| Control plants | 3,749 |
| Case median predicted pChEMBL | 4.876 |
| Control median predicted pChEMBL | 4.886 |
| Median shift (Hodges-Lehmann) | -0.003 |
| Family-clustered 95% CI | -0.027 to +0.035 |
| Mann-Whitney p | 0.8373 |
| Rank-biserial correlation | -0.011 |

## Negative control

The same comparison using a traditional-use category with no expected PTP1B relevance. If this matched the primary effect, the primary result would be attributable to study bias rather than pharmacological convergence.

- Negative-control plants: 54
- Median shift: +0.013
- Mann-Whitney p: 0.4574

## Prespecified support rule

| Condition | Met |
|---|:--:|
| p below alpha | no |
| ci excludes zero | no |
| negative control smaller | no |

All three conditions were required. Outcome: **NOT SUPPORTED**.

## Cohort

- Duke ethnobotanical use records: 82,873
- NPASS source organisms: 48,940
- Plants covered by both sources: 4,388
- Plants analysed: 3,870 across 305 families
- Compounds scored: 111,485
- Compounds unusable (absent or unparsable structure): 654

## Limitations

- Predicted activity is not measured activity. This tests convergence with a model, not experimental confirmation.
- Duke coverage is biased toward well-studied and temperate/North American floras.
- Traditional 'antidiabetic' use is recorded against a modern disease label applied retrospectively to heterogeneous historical sources.
- NPASS compound-organism edges are reported provenance, not independently verified occurrence.
- A single target (PTP1B) cannot represent antidiabetic activity in general.
- Exact species-name matching misses synonyms and will undercount true overlap.
- This analysis does not and cannot establish efficacy or safety of any plant or preparation.

## Traditional knowledge handling

- Every traditional-use claim retains its Duke REFERENCE and COUNTRY fields in the output tables.
- Convergence evidence only. The analysis must never state that traditional knowledge is validated or invalidated by a model prediction.
- No genetic resources are accessed. Only published literature records are used.

## Input checksums

| Input | SHA-256 |
|---|---|
| species_info | `d7d4d429e64de1268ba9d5c69f7d72fe16674731b04362b239f801b63d17d275` |
| species_pair | `0c80055682df5e648941833876c03956cbb4fa3254d1ac6d5dd1f2896bd480cb` |
| structures | `5f2ca5c3fcd0386d65a31146122e4c8bf6dfcbea7a786ea01e6be764d69f65dc` |
| duke | `ff325bc477c355799d12194e2d172c9e181ba44e59056c5ca97f41db595b7d19` |

