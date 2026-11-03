# PhytoTRUST study-connected ChEMBL dataset

## Construction contract

The candidate list and identifier-only selection criteria were frozen before download. 12510 uncapped ChEMBL IC50 records were retrieved from 10 previously unused targets. Previously locked evaluation outcomes were never used to rank targets or contexts; records sharing a locked compound, activity, assay, or publication identifier were excluded.

Only exact, positive human IC50 values for ChEMBL SINGLE PROTEIN targets were retained. Reported and computed pChEMBL had to agree within 0.15, data-validity comments and potential duplicates were excluded, structures had to standardize successfully, and endpoint/organism/assay-type/BAO context was never pooled across context keys.

## Frozen target selection

| Target | Eligible rows | Compounds | Studies | Repeated compounds | Multi-compound studies | Retained |
|---|---:|---:|---:|---:|---:|---:|
| CHEMBL222 | 1602 | 1414 | 96 | 53 | 89 | True |
| CHEMBL233 | 204 | 188 | 8 | 14 | 7 | False |
| CHEMBL218 | 94 | 82 | 10 | 10 | 7 | False |
| CHEMBL236 | 310 | 300 | 38 | 8 | 31 | False |
| CHEMBL214 | 136 | 119 | 18 | 7 | 15 | False |
| CHEMBL217 | 313 | 241 | 31 | 6 | 24 | False |
| CHEMBL219 | 109 | 99 | 22 | 2 | 11 | False |
| CHEMBL224 | 223 | 208 | 59 | 2 | 27 | False |
| CHEMBL225 | 129 | 126 | 29 | 2 | 16 | False |
| CHEMBL226 | 87 | 82 | 18 | 1 | 7 | False |

## Final graph and partitions

| Target | Partition | Measurements | Edges | Compounds | Studies | Components | Largest component C/S | Repeated compounds | Studies >=3 compounds | 4-cycles | Verified NP compounds | Structural gate |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CHEMBL222 | development | 966 | 907 | 877 | 61 | 49 | 88/5 | 29 | 52 | 78 | 10 | True |
| CHEMBL222 | future_evaluation | 636 | 591 | 537 | 35 | 4 | 500/32 | 24 | 34 | 69 | 21 | True |

## Leakage and provenance

Development and future evaluation are disjoint by connected compound-publication component. The final dataset contains 1602 measurements, 1414 compounds, 96 publications, 1498 unique compound-publication edges, and 31 structure-verified NPASS/COCONUT compounds (27 COCONUT, 13 NPASS, 9 in both). Compound, publication, assay, and activity overlap between partitions are all zero.

Future-evaluation components contain at least one publication from 2020 onward. Their pre-2020 rows are saved separately as permissible historical evidence; 2020+ values are in the explicitly quarantined labels file and were not summarized or used in connectivity selection.

Every measurement retains ChEMBL activity, assay, publication, molecule, structure, context, original units, pChEMBL, snapshot checksum, retrieval time, and direct ChEMBL URLs. NP membership retains database source, source IDs, versions, and standardized structural identity.

## Identifiability verdict

This graph gate is necessary but not sufficient for statistical identifiability; outcome-based variance estimation is intentionally outside this phase.

Development targets passing the necessary structural gate: CHEMBL222. Targets failing it: none.

The dataset supports a limited single-target random-intercept feasibility study, not a multi-target or journal-level generalization claim. Only one frozen candidate passed selection; its development graph has 29 compounds repeated across studies, only one compound in at least three studies, and 49 disconnected components.

Passing this gate does not prove that a variance component will be estimable or predictive. It only establishes repeated compounds, multi-compound studies, and graph cycles needed to separate some compound and study effects. Outcome-based variance estimation is intentionally deferred to a later modeling phase.

## Limitations

- ChEMBL publication identifiers are not guaranteed to represent independent laboratories; laboratory identity is usually unavailable.
- BAO format and assay descriptions do not capture substrate concentration, construct, temperature, pH, incubation time, or all biological conditions. Missing conditions remain unknown.
- The selected connected contexts may not represent the full target literature or real-world disagreement prevalence.
- Natural-product status is exact NPASS/COCONUT structural membership, not proof of the provenance of the assayed sample, efficacy, or safety.
- Component-disjoint partitioning is conservative and can place substantial pre-2020 history into the future-evaluation component; those rows must not be used to fit development models except as per-case historical input under a frozen evaluation protocol.
- No M0/M1 model was trained, fitted, or evaluated in this phase.
