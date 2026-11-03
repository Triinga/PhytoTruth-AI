# Evidence Reversal Experiment

## Scientific question

Does a prespecified evidence-aware ranking prioritize compounds whose experimental potency reproduces in later publications better than ranking by the strongest earlier measurement?

This is a retrospective evidence-replication benchmark. It is not a clinical or prospective validation, and its reliability score is a transparent heuristic rather than a calibrated probability.

## Prespecified cohort

- Target: ChEMBL `CHEMBL335`.
- Endpoint: exact IC50 values normalized to pActivity.
- Organism: *Homo sapiens*.
- Assay type: binding (`B`).
- BAO format: single-protein format (`BAO_0000357`).
- Compatible context selection: the context with the largest number of distinct eligible molecules, without using activity outcomes.
- Natural-product status: standardized-structure match to the complete ingested NPASS 3.0 structure release. Name-based inference is forbidden.

Censored measurements, missing publication years, duplicate activity IDs, non-finite activities, and other assay contexts are excluded from the primary benchmark and counted in `results.json`.

## Leakage-resistant splits

### Primary: latest-publication-year holdout

For each molecule with measurements in at least two publication years, all records from its latest year are held out. Only strictly earlier records may contribute to any ranking score. Activity IDs and document IDs must be disjoint. This yields more independent molecules than a single global cutoff while retaining strict temporal ordering within every molecule.

### Sensitivity: global publication cutoff

The prespecified cutoff is 2015. Records published through 2015 score the molecule; records after 2015 are held out. Only molecules observed on both sides are evaluated. This more closely approximates a historical prospective test but has a smaller cohort.

The unit for uncertainty and statistical inference is the molecule, not the individual measurement.

## Rankings

- Baseline: maximum pActivity among training records.
- Median potency: median training pActivity.
- PhytoTruth full: the existing evidence score, with explicit penalties for heuristic unreliability and bootstrap-interval width.
- Ablations: remove reliability, uncertainty, explicit consistency, or support one at a time. No weight is fitted or selected using held-out outcomes.

The full score remains a ranking heuristic. Values must not be interpreted as probabilities of efficacy, reproducibility, safety, or clinical success.

## Held-out outcome and metrics

The outcome is median pActivity in held-out evidence. The activity threshold is pActivity 6.0 and K is 10. Reported metrics are nDCG@10 with linear log-potency gain above the threshold, precision@10, enrichment factor@10, and Kendall tau-b.

Uncertainty uses 2,000 paired molecule-level bootstrap samples. Method comparisons use 2,000 paired score-swap randomizations with deterministic seed 2026. Confidence intervals and p-values are exploratory because several related endpoints and cohorts are reported; there is no multiplicity-adjusted confirmatory claim.

## Reproduction

After running the documented data pipeline:

```powershell
cd backend
python -m app.pipeline evidence-reversal `
  --activities ../data/processed/chembl335_ic50.normalized.jsonl `
  --natural-products ../data/processed/chembl335_ic50.natural-products.jsonl `
  --output-dir ../experiments/evidence_reversal_chembl335 `
  --global-cutoff-year 2015 --k 10 --active-threshold 6.0 `
  --bootstrap-iterations 2000 --permutation-iterations 2000 --seed 2026
```

The output directory contains input checksums, the exact command, cohort CSV files, machine-readable results, a metric summary, an SVG comparison plot, and an automatically generated report.

## Executed result (2026-10-20)

- Compatible input: 2,915 measurements from 2,451 molecules.
- Latest-year cohort: 162 molecules; NPASS-verified subgroup: 31.
- Global-2015 cohort: 39 molecules; NPASS-verified subgroup: 14.
- All-compound latest-year nDCG@10: baseline 0.891, PhytoTruth 0.977; paired difference 0.086, 95% CI 0.000 to 0.265, p=0.3458.
- Natural-product latest-year nDCG@10: baseline 0.891, PhytoTruth 0.765; difference -0.126, 95% CI -0.758 to 0.418.
- Median potency alone matched the full method's top-10 metrics and had higher Kendall tau-b. Full-minus-median Kendall difference was -0.035, 95% CI -0.079 to -0.003, p=0.0400.

Conclusion: this experiment does **not** demonstrate that the current evidence-aware heuristic outperforms potency-only ranking. The all-compound top-10 point estimate improved without statistically supported paired evidence, the natural-product point estimate worsened, and the ablation suggests the benefit came from replacing maximum potency with robust median potency rather than from the reliability components. This negative result should guide future redesign and prospective validation; it must not be presented as a successful validation.

## Limitations

- Only 31 NPASS-verified natural products met the repeated, temporally separated evidence criteria in the primary compatible context.
- The global temporal cohort is smaller still, so intervals are wide.
- Metadata compatibility does not guarantee identical laboratory protocols.
- Literature measurements can share unobserved dependencies even when document IDs differ.
- The benchmark covers one target and one endpoint; results may not generalize.
- COCONUT was not part of the complete membership index used for this run, so natural-product verification is NPASS-only.
