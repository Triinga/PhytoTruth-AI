# Independent assay-disagreement validation

## Scientific question

Can the already-frozen assay-disagreement models predict a ten-fold or larger shift in exact IC50 potency for compounds, targets, assays, and publications that were never used in earlier development or evaluation?

This phase is evaluation only. It does not retrain, recalibrate, tune, select, or modify a model.

## Frozen design

`experiments/independent_validation_protocol.json` was written before disagreement outcomes were constructed. Five previously unused human single-protein targets were selected from biological relevance and ChEMBL exact-IC50 volume: VEGFR2, estrogen receptor, androgen receptor, HER2, and ABL1. The original model targets were acetylcholinesterase, EGFR, COX-2, DPP-4, and PTP1B.

The validation endpoint is unchanged: the absolute difference between the median pActivity of all earlier-publication exact-IC50 evidence and the median pActivity of the latest publication year, with disagreement at least 1.0 log10 unit. Earlier and held-out assays and publications must be disjoint. Held-out years must be 2020 or later.

Prespecified overall minimums are 100 events, 20 disagreements, and 10 held-out publications. Natural-product minimums are 30 events, 10 disagreements, and five publications. The support rule requires the frozen full model's study-clustered AUROC lower bound above 0.5, PR-AUC lower bound above observed prevalence, and Brier score no worse than the frozen constant-prevalence fallback.

## Source retrieval and preparation

Use the official ChEMBL REST API and the existing NPASS/COCONUT index. No record cap is used.

```powershell
cd backend
$targets = @('CHEMBL279','CHEMBL206','CHEMBL1871','CHEMBL1824','CHEMBL1862')
foreach ($target in $targets) {
  python -m app.pipeline fetch-chembl --target $target --standard-type IC50 `
    --output "../data/independent_validation/$target.raw.jsonl" --page-size 1000
  python -m app.pipeline normalize-chembl `
    --input "../data/independent_validation/$target.raw.jsonl" `
    --output "../data/independent_validation/$target.normalized.jsonl" `
    --rejects "../data/independent_validation/$target.rejects.jsonl"
  python -m app.pipeline crossmatch-natural-products `
    --input "../data/independent_validation/$target.normalized.jsonl" `
    --database "../data/processed/natural_products.sqlite" `
    --output "../data/independent_validation/$target.natural.jsonl" `
    --rejects "../data/independent_validation/$target.natural.rejects.jsonl"
}
```

Raw, normalized, and natural-product files have adjacent checksum-bearing manifests. NPASS/COCONUT membership is determined by standardized structure and is used only for subgroup reporting.

## Leakage controls

The exclusion ledger reads every record—not just modeled events—from all previous five-target normalized snapshots. It contains 27,996 compounds, 41,274 activities, 4,415 assays, and 2,736 publications. A new event is rejected if any compound, activity, assay, publication, or target overlaps that ledger. The previous 531-event dataset is indexed again as defense in depth.

Every compound must occur once in the final validation table. If a compound creates events for multiple new targets, all such events are excluded rather than choosing a favorable target. Automated assertions verify the temporal boundary and all identifier intersections before scoring.

## Frozen models

- `assay_disagreement_model.joblib`: the original full logistic model.
- `study_generalized_model.joblib`: the development-selected constant-prevalence fallback.

Both files are SHA-256 hashed before and after prediction. Scoring aborts if a hash changes. The five validation targets are unknown to the fitted one-hot encoder, so this deliberately measures cross-target and cross-study transfer.

## Inference and outputs

AUROC, PR-AUC, Brier score, log loss, ECE-10, calibration intercept, and calibration slope are calculated without refitting. Confidence intervals resample held-out publication clusters within target. Natural products are evaluated separately only when both outcome classes exist.

Outputs include:

- `independent_dataset.csv`: one row per accepted compound event;
- `validation_measurements.csv`: every real prior/held-out measurement used, linked to its event and role;
- `frozen_model_predictions.csv`: untouched model scores and provenance;
- `independence_exclusions.csv`: auditable rejected candidates and reasons;
- `target_metrics.csv`, `results.json`, an SVG comparison, report, command, and checksum manifest.

## Executed result and limitations

The frozen target set produced 263 candidate repeated-measurement events. Strict exclusions left 90 events from 37 held-out publications with seven disagreements. Only five were verified natural products and none disagreed, so NP AUROC and PR-AUC are undefined.

The planned overall minimum was missed by 10 events and 13 disagreements; the natural-product minimum was missed by 25 events and 10 disagreements. These are count deficits, not a power calculation. Adding records specifically because they disagree would bias prevalence; additional data must come from prospectively frozen targets or future ChEMBL releases before outcomes are examined.

The full model had AUROC 0.494, PR-AUC 0.093, Brier 0.421, ECE-10 0.538, calibration intercept -2.448, and calibration slope -0.068. The frozen fallback had Brier 0.097. The evidence does not support reliable independent prediction.
