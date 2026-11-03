# Assay protocol metadata enrichment report

## Scope and integrity controls

This release enriches **531** assay-disagreement events with **1050** exact ChEMBL assay records. The raw API response, request URL, retrieval timestamp, and record SHA-256 are retained for every assay. Extraction version: `assay-protocol-rules-v1`.

Structured ChEMBL fields are copied without reinterpretation. Description-derived categories require literal rule matches and include the matched pattern and source-text snippet. Empty categories mean *not reported or not detected*, never that a method or condition was experimentally absent.

Prior-assay fields are labelled `pre_prediction_available`. Held-out protocol fields and prior-versus-held-out comparisons are labelled `post_outcome_explanatory_only`; they must not be used to predict the held-out result.

All included assays are ChEMBL Binding assays in single-protein BAO format, reflecting the benchmark's compatibility filter. Those two fields have no within-cohort variation and therefore cannot explain disagreement here.

## Metadata coverage

| Field | Assays reported / explicitly mentioned | Coverage |
|---|---:|---:|
| assay_type_description | 1050 / 1050 | 100.0% |
| assay_organism | 807 / 1050 | 76.9% |
| assay_tissue | 12 / 1050 | 1.1% |
| assay_cell_type | 0 / 1050 | 0.0% |
| assay_subcellular_fraction | 0 / 1050 | 0.0% |
| assay_strain | 0 / 1050 | 0.0% |
| bao_label | 1050 / 1050 | 100.0% |
| description | 1050 / 1050 | 100.0% |
| assay_category | 5 / 1050 | 0.5% |
| assay_test_type | 5 / 1050 | 0.5% |
| assay_parameter_types | 13 / 1050 | 1.2% |
| variant_mutation | 118 / 1050 | 11.2% |
| methodology_categories | 132 / 1050 | 12.6% |
| detection_categories | 257 / 1050 | 24.5% |
| biological_system_text_categories | 221 / 1050 | 21.0% |
| endpoint_categories | 123 / 1050 | 11.7% |
| temperature_mentions | 29 / 1050 | 2.8% |
| ph_mentions | 44 / 1050 | 4.2% |
| duration_mentions | 248 / 1050 | 23.6% |
| substrate_mentions_json | 358 / 1050 | 34.1% |

## Dataset partitions

| Partition | Events | Disagreement rate | Verified natural products | Complete prior descriptions | Complete held-out descriptions |
|---|---:|---:|---:|---:|---:|
| development | 212 | 23.6% | 47 | 212 | 212 |
| evaluation | 319 | 53.0% | 60 | 319 | 319 |

## Protocol differences and disagreement

Risk differences below compare observed disagreement rates for events whose prior and held-out protocol metadata differed versus matched. Confidence intervals use held-out-study cluster bootstrap. Fisher tests and BH q-values are supplementary because molecule rows within a study are not independent. These are retrospective associations, not causal effects and not model validation.

| Evaluation feature | Comparable | Changed | Risk difference | Study-clustered 95% CI | BH q |
|---|---:|---:|---:|---:|---:|
| protocol_biological_system_text_changed | 19 | 4 | +0.483 | [-0.353, +0.882] | 0.534 |
| protocol_duration_changed | 29 | 26 | -0.359 | [-0.815, +0.296] | 0.534 |
| protocol_ph_changed | 18 | 7 | -0.338 | [-0.338, -0.338] | 0.534 |
| protocol_organism_changed | 290 | 3 | +0.137 | [-0.560, +0.710] | 1 |
| protocol_detection_changed | 32 | 23 | +0.082 | [-0.274, +0.445] | 1 |
| protocol_variant_mutation_changed | 31 | 24 | +0.048 | [-0.441, +0.600] | 1 |

### Prediction-time metadata missingness

These fixed comparisons use only prior-assay metadata available before the held-out study. A positive risk difference means events with reported prior metadata had more later disagreement than events where that field was unreported; it does not imply the metadata caused disagreement.

| Evaluation feature | Comparable | Reported | Risk difference (reported - missing) | Study-clustered 95% CI | BH q |
|---|---:|---:|---:|---:|---:|
| prior_detection_reported | 319 | 80 | -0.290 | [-0.459, +0.052] | 9.64e-05 |
| prior_biological_system_text_reported | 319 | 65 | -0.279 | [-0.443, +0.047] | 0.000439 |
| prior_endpoint_reported | 319 | 53 | -0.273 | [-0.451, +0.063] | 0.00105 |
| prior_variant_reported | 319 | 46 | -0.263 | [-0.449, +0.095] | 0.00341 |
| prior_tissue_reported | 319 | 12 | -0.204 | [-0.451, +0.210] | 0.438 |
| prior_parameters_reported | 319 | 3 | -0.198 | [-0.638, +0.618] | 0.737 |
| prior_duration_reported | 319 | 67 | -0.198 | [-0.384, +0.125] | 0.0124 |
| prior_ph_reported | 319 | 35 | +0.079 | [-0.077, +0.307] | 0.651 |
| prior_methodology_reported | 319 | 65 | -0.066 | [-0.167, +0.150] | 0.634 |
| prior_organism_reported | 319 | 292 | -0.028 | [-0.411, +0.302] | 0.926 |
| prior_temperature_reported | 319 | 24 | +0.013 | [-0.304, +0.188] | 1 |

## Interpretation

No protocol-difference or prior-metadata-availability association met a conservative transfer criterion of directional agreement in development, an evaluation study-clustered interval excluding zero, and adequate independent-study support. Several unclustered Fisher/BH results appear small, but their study-clustered intervals cross zero; they must not be treated as validated signals.

The defensible finding is therefore about data adequacy: ChEMBL reliably supplies an assay description, broad assay type, and BAO format for this cohort, but most experimental conditions are too sparsely and unevenly reported to establish that protocol differences explain disagreement. Full-text methods and more independent natural-product studies are required.

Natural-product subgroup rows are retained in `protocol_associations.csv`; sparse comparable groups and study concentration must be considered before interpretation.

## Limitations

- ChEMBL assay records frequently omit temperature, pH, duration, substrate, detection method, cell type, tissue, or strain. Free-text extraction cannot recover unreported conditions.
- Literal keyword categories improve auditability but do not replace manual protocol curation or full-text methods sections.
- ChEMBL descriptions can summarize rather than reproduce the experimental protocol. The raw record is authoritative for this artifact, not the derived category.
- Protocol comparisons use metadata available only after the held-out assay was performed and are explanatory only. No prediction model was trained or retrained in this phase.
- Associations can be confounded by target, era, publication practice, and study composition; they do not establish that a protocol difference caused disagreement.

## Reproduction

```powershell
python -m app.pipeline fetch-assay-metadata --events ../experiments/assay_disagreement_generalization/dataset.csv --output ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl
python -m app.pipeline enrich-assay-protocols --events ../experiments/assay_disagreement_generalization/dataset.csv --assays ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl --output-dir ../experiments/assay_protocol_enrichment --bootstrap-iterations 1000 --seed 2026
```

Input assay archive SHA-256: `97262a32685a812851225f241a438f96df3fc99f67b2607f2de65805ecc872b9`.
