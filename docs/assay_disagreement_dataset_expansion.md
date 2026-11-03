# Assay Disagreement Dataset Expansion

## Scope and frozen design

This phase expands the independent assay-disagreement data only. It does not
load, retrain, tune, or evaluate any predictive model. Fifteen human
single-protein targets were prespecified before disagreement outcomes were
constructed. Non-human homologues were not pooled with human targets.

An event is one compound-target pair with exact IC50 measurements in at least
two publication years. Historical evidence comprises all eligible earlier years;
the independent evidence is the latest year. The two sides must have disjoint
ChEMBL document and assay identifiers. A positive label requires an absolute
difference of at least 1.0 between the historical and independent median
pChEMBL values. This operational threshold was frozen before outcome review.

Within each target, a single context (target, target organism, endpoint, broad
assay type, and BAO format) is selected by eligible molecule count without using
activity values. Exact relations are required. Censored records, ChEMBL
potential duplicates, validity-comment records, inconsistent pChEMBL
normalization, incompatible endpoints, and other contexts do not become labels.

## Data obtained

- ChEMBL API snapshot: 49,632 uncapped IC50 records across 15 targets.
- Normalization: 49,632 accepted, 0 normalization rejects; 60 records lacked a
  structure for natural-product matching.
- Exact NPASS/COCONUT structural matches: 4,845 measurement rows before event
  construction.
- Repeated-study candidates before strict prior-exposure filtering: 394.
- New accepted events: 94 unique compounds, 58 held-out studies, 12 targets.
- Labels: 24 positive disagreements and 70 negatives.
- Natural products: 20 events and 7 positive disagreements. Thirteen compounds
  matched COCONUT only, six matched both COCONUT and NPASS, and one matched
  NPASS only.
- Used evidence: 295 ChEMBL measurements and 196 distinct assays; all 196 assay
  metadata records were returned by the ChEMBL API.

Together with the earlier 90-event independent validation, the project now has
184 non-overlapping independent events, 31 positives, and 25 verified
natural-product events. The expansion remains a separate, quarantined dataset;
the combined count is descriptive and no performance evaluation was performed.

## Leakage and exclusions

The exposure union contains 10 earlier targets, 44,750 molecules, 71,418
activities, 7,217 assays, and 4,009 documents. An accepted event may overlap none
of them. It must also be from 2020 or later and have a molecule not repeated in
another new target event.

The strict candidate filter rejected 300 of 394 candidate events. Reasons are
not mutually exclusive: 208 were before the frozen temporal boundary (and hence
not in the evaluation partition), 110 reused an earlier publication, 88 reused
an earlier compound, and 58 involved compounds repeated across new targets.
There was no earlier activity, earlier assay, or within-event assay/document
overlap among accepted records. Before candidate construction, 4,508 records
were excluded because ChEMBL flagged them as potential duplicates and 23,395
single-publication-year compounds could not form independent-study events.

## Protocol interpretation

All 94 accepted events agree on broad assay type and BAO format by construction.
Assay organism agrees for 91 and is unreported on one or both sides for 3; no
accepted event uses a documented organism change as a disagreement. More
detailed fields are sparse: parameter type, pH, strain, subcellular fraction,
temperature, tissue, and variant comparisons are unknown for every event.

Among 24 positive cases:

- 4 have a documented protocol difference that could explain the shift: one
  detection-method difference, one endpoint-description difference, and two
  duration differences.
- 3 have at least four reported dimensions that agree, but remain only partially
  documented and therefore are not proof of unexplained biological variation.
- 17 cannot be interpreted because protocol metadata are incomplete.

These are descriptive post-outcome annotations, not predictive features and not
causal findings. The JSON case-study artifact lists every historical and
independent activity, value, assay, document, year, structure, natural-product
source, shift magnitude, documented comparison, and the remaining uncertainty.

## Reproduction

The exact commands are in the root README. Raw ChEMBL activity and assay
snapshots have retrieval timestamps and SHA-256 manifests. Generated output
checksums are recorded in `experiments/assay_disagreement_expansion/final/manifest.json`.
The integrity tests enforce the frozen threshold, exact IC50 relations, event ID
determinism, unique molecules, temporal boundary, assay/document disjointness,
complete measurement linkage, natural-product provenance, and label quarantine.

## Limitations

- ChEMBL assay descriptions often summarize rather than reproduce full methods.
- The largest compatible context is deliberately outcome-blind but may omit
  informative smaller assay formats.
- A latest-year event may contain several activity rows from one publication;
  those rows are not treated as independent studies.
- NPASS/COCONUT matching verifies standardized structural membership only.
- Twenty natural products and seven natural-product positives remain too few for
  reliable subgroup model evaluation; further independent publications or
  manually curated full-text assay methods are required.
