# Assay protocol metadata enrichment

## Purpose

This phase asks whether omitted protocol metadata could help explain the assay-disagreement model's poor cross-study generalization. It enriches the already frozen disagreement events; it does not change their measurements, partitions, outcomes, natural-product labels, or model results.

## Data source and retrieval

Assay identifiers are the union of `prior_assay_ids` and `heldout_assay_ids` in `experiments/assay_disagreement_generalization/dataset.csv`. `fetch-assay-metadata` sorts and deduplicates those identifiers and queries the official ChEMBL REST `assay` endpoint in bounded batches. HTTP 408, 425, 429, and 5xx responses use the existing exponential-backoff retry mechanism.

Each JSONL envelope contains:

- the complete ChEMBL assay response without field removal;
- the exact request URL and UTC retrieval timestamp;
- the assay ChEMBL identifier; and
- a SHA-256 of the canonicalized raw assay object.

The adjacent manifest records the frozen event-dataset checksum, requested and returned counts, missing identifiers, fetch settings, and archive checksum. Retrieval is atomic through a `.part` file.

## Structured fields

The feature table copies ChEMBL's controlled or structured values directly, including assay type, target, document, organism/taxonomy, tissue, cell type, subcellular fraction, strain, BAO format, assay category/test type/group, confidence and relationship metadata, assay parameters, and variant accession/mutation/organism. `null`, an empty parameter list, or absent variant object remains missing.

The full raw parameter and variant records remain in the raw JSONL. The CSV carries compact exact parameter types and variant identifiers instead of duplicating full protein sequences.

## Description extraction

`assay-protocol-rules-v1` is deterministic and case-insensitive. It emits a category only after a literal regular-expression match in the ChEMBL `description` field. Fixed families cover:

- methodology: binding, enzyme activity, functional, cell-based, viability/proliferation, phosphorylation;
- detection: fluorescence, luminescence, radiometric, absorbance/colorimetric, FRET, AlphaScreen, scintillation proximity, immunoassay, chromatography/mass spectrometry;
- biological-system text: recombinant, purified, cellular, membrane, lysate, tissue;
- endpoints: IC50, EC50, GI50, Ki, Kd, percent inhibition, activity; and
- explicit temperature, pH, incubation/reaction duration, and `substrate` text mentions.

Every derived row retains the rule version, source field, matched pattern, and verbatim surrounding snippet as JSON. The extractor does not infer a substrate from target identity, a detection method from assay type, a species from a protein accession, or an endpoint from the downstream benchmark filter. A non-match is therefore `not reported or not detected`, never `absent`.

## Timing and leakage controls

- `prior_*` protocol summaries are labelled `pre_prediction_available` because they refer only to evidence published before the held-out experiment.
- `heldout_*` fields are labelled `post_outcome_explanatory_only`.
- `protocol_*_changed` flags compare prior with held-out protocols and are also `post_outcome_explanatory_only`.

The comparison flags cannot be used as predictors of the assay they describe. This phase performs no model selection, fitting, calibration, or threshold tuning.

## Descriptive analysis

Analysis is fixed before looking at associations. For each development/evaluation partition and for all/verified-natural-product subgroups, it computes the observed disagreement risk difference between known-changed and known-same protocol sets. Unknown comparisons are excluded rather than treated as same. Confidence intervals resample held-out publication clusters, preserving within-study dependence. Fisher exact p-values with Benjamini-Hochberg correction are supplementary row-level summaries; clustered intervals are primary because compounds from a publication are correlated.

A protocol association is considered descriptive, not causal. A scientifically persuasive transfer finding would require adequate comparable events, consistent direction in development and untouched evaluation, a study-clustered interval excluding zero, and replication across targets/studies. Sparse natural-product rows do not satisfy those conditions.

## Outputs

- `assay_metadata.raw.jsonl` and manifest: immutable source archive and retrieval provenance;
- `assay_protocol_features.csv`: one row per ChEMBL assay;
- `enriched_assay_disagreement_events.csv`: original events plus timing-labelled summaries;
- `coverage_statistics.json`: field and partition coverage;
- `protocol_associations.csv` / `.json`: fixed descriptive comparisons;
- `prior_metadata_missingness_associations.csv` / `.json`: prediction-time metadata availability versus later disagreement;
- `heldout_study_metadata_coverage.csv`: study-level coverage and missingness patterns;
- `protocol_coverage_and_associations.svg`: deterministic summary figure;
- `SCIENTIFIC_REPORT.md`: automatically generated findings and limitations; and
- `manifest.json`: input/output hashes, parameters, seed, and timing policy.

## Limitations

ChEMBL descriptions are summaries and rarely include a complete methods protocol. Full-text articles or supplementary methods, curated with document-level licensing and provenance, would be required for systematic buffer composition, reagent lots, exact substrate concentrations, plate formats, incubation sequences, and instrument settings. More independent natural-product studies are needed before subgroup conclusions are reliable.
