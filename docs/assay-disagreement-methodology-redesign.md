# Assay disagreement methodology redesign

This phase asks whether expected assay-to-assay variation is more learnable than a thresholded disagreement label. It is deliberately confined to the 212 events dated through 2019 in the original development partition. The prior independent cohorts (90 and 94 cases) remain locked and are neither loaded nor evaluated.

## Design

- Primary outcome: absolute change between the median prior and held-out exact-IC50 pChEMBL measurements.
- Secondary outcome: at least 1.0 pChEMBL absolute change.
- Comparators: fold-specific constant prevalence/median, a Jeffreys-smoothed historical contradiction score, and the prior logistic architecture.
- Redesigned methods: robust Huber regression and random-forest magnitude regression. Tree votes above 1.0 provide a risk score; tree quantiles provide a descriptive ensemble interval, not a calibrated confidence interval.
- Feature ablation: historical evidence alone versus historical evidence plus molecular descriptors and prior-only protocol-reporting indicators.
- Splits: connected publication–compound groups, forward publication time, and complete target holdout. Natural-product status is never a feature.

The exact seed, methods, numerical settings, amendment, validation rules, and readiness rule are frozen in `experiments/methodology_redesign_protocol.json`. The post-first-run amendment only corrects Huber non-convergence and unconstrained extrapolation; it does not change the primary candidate or readiness rule.

## Reproduction

```powershell
cd backend
python -m app.methodology_redesign --dataset ../experiments/assay_protocol_enrichment/enriched_assay_disagreement_events.csv --output-dir ../experiments/assay_disagreement_methodology_redesign --bootstrap-iterations 1000 --seed 2026
python -m pytest tests/test_methodology_redesign.py -q
```

The generated `SCIENTIFIC_REPORT.md` is authoritative for executed counts and results. `development_predictions.csv` retains ChEMBL document and assay provenance; `protocol_ablation.csv` contains paired study-cluster bootstrap comparisons; `natural_product_error_analysis.csv` exposes compound-level failures.

## Interpretation limits

The risk scores are research scores, not probabilities of biological failure or clinical inefficacy. Protocol missingness indicators do not recover unreported experimental conditions. No candidate from this phase is ready for deployment or a new independent evaluation unless every prespecified development-readiness condition is met.
