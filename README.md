# PhytoTruth AI

**Evidence before excitement.** PhytoTruth is an end-to-end research MVP for EthnoHack 2026 Track 2. It combines real natural-product source adapters, experimental bioactivity evidence, RDKit fingerprints, leakage-resistant ML validation and a modern evidence-audit interface.

> Research prioritization only. The system does not make medical recommendations or claims of efficacy or safety.

## What is working

- A FastAPI service exposes a real, traceable ChEMBL PTP1B demonstration slice retrieved on 2026-10-20, including activity, assay and document IDs.
- The evidence engine detects potency reversals, large assay spread, duplicates, source validity comments and organism mixing; it reports bootstrap uncertainty.
- RDKit generates canonical structures, Morgan radius-2 fingerprints, descriptors, Murcko scaffolds and 2D molecule SVGs.
- Random Forest training uses scaffold-grouped cross-validation.
- `POST /api/predict` serves the selected real XGBoost model from checksum-validated native artifacts, with retrospective uncertainty and an explicit applicability-domain warning.
- The robust ranking experiment compares development-selected reliability-aware ranking with median potency across whole held-out biological targets and future-publication evidence.
- Streaming COCONUT and NPASS CSV/TSV/SDF importers validate and standardize structures with RDKit, deduplicate by standard InChIKey, and retain release, source ID, file, line, original SMILES, and raw-row provenance in a disk-backed SQLite index.
- React + TypeScript + Tailwind provides responsive discovery, evidence, ranking-experiment and data-lineage workspaces.
- The prediction workspace calls the real `/api/predict` model, renders RDKit structures, descriptors, residual intervals, and explicit applicability-domain warnings without substituting demo values.
- A semantic SQLite knowledge graph connects reported source organisms, structure-verified compounds, ChEMBL targets, curated diseases, assays and publications. Every edge carries an evidence source and identifier; unsupported biological or clinical relationships are never inferred.
- A bounded-memory, multi-target screening pipeline streams the complete NPASS/COCONUT index, applies configurable physicochemical filters and PAINS/Brenk research alerts, scores checksum-validated target models with applicability-domain and uncertainty penalties, and exports reproducible ranked candidates to CSV/JSONL plus the graph.
- The Reversal Lab reads the generated multi-target artifact and exposes target-cluster confidence intervals, exact target sign-flip p-values, effect sizes, power limits, and ablations—including the negative primary conclusion.
- The Study Planner converts a robust Student-t future-study distribution into disagreement risk, uncertainty-based abstention, and outcome-blind next-experiment recommendations based on expected information gain.

The bundled five-molecule dashboard is deliberately a small evidence demonstration, not a claim that all five molecules are natural products. Okadaic acid and calyculin A are labeled as known natural products; unlabeled compounds remain `not yet cross-matched` until the COCONUT/NPASS importer establishes source provenance.

## Repository

```text
backend/app/          FastAPI, evidence scoring, RDKit and ML pipeline
backend/data/         checked-in real-source demo snapshot
backend/tests/        reliability-engine tests
experiments/          reproducible experiment and inference artifacts
frontend/src/         React research workspace
data/raw/             ignored immutable source downloads
data/processed/       ignored normalized artifacts
docs/                 architecture and experimental protocol
models/               ignored trained model artifacts
```

## Run locally

Recommended: Python 3.11 or 3.12 and Node 20+.

```powershell
cd backend
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
uvicorn app.main:app --reload
```

In another terminal:

```powershell
cd frontend
npm install
npm run dev
```

Open `http://localhost:5173`. If the API is absent, the interface uses a clearly indicated checked-in ChEMBL snapshot; molecule rendering requires the API and RDKit.

The API's interactive documentation is at `http://127.0.0.1:8000/docs`. A prediction request is:

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/api/predict `
  -ContentType application/json `
  -Body '{"smiles":"CCO"}'
```

Inference uses `experiments/ml_activity_chembl335` by default. Override this with `PHYTOTRUTH_MODEL_DIR`. See the [prediction API guide](docs/api.md) for artifact conversion, complete response fields, errors, and scientific limitations.

### Frontend scientific dashboard

The five workspaces deliberately separate different evidence classes:

- **Discovery** and **Evidence** show measured ChEMBL records and deterministic reliability flags.
- **Predict** sends the submitted SMILES to the saved XGBoost model. Results are labeled `ML PREDICTION`; out-of-domain structures are marked unreliable.
- **Reversal Lab** visualizes `experiments/robust_multitarget_evaluation/results.json`. It keeps median potency deployed because the optimized ranker did not demonstrate superiority.
- **Study Planner** visualizes `experiments/evidence_intelligence/results.json`. Its experiment recommendations are research-planning hypotheses and never use held-out outcomes for priority.
- **Data Sources** reports the actual integration state, including the complete NPASS structure index and the bounded COCONUT parser verification.

Run frontend checks and create a production bundle:

```powershell
cd frontend
npm install
npm test
npm run build
```

The `predev`, `pretest`, and `prebuild` hooks copy the executed single-target, robust multi-target, and study-planning result artifacts into the frontend bundle and fail if an artifact is absent. This prevents the dashboard from silently substituting invented or stale metrics.

Set `VITE_API_URL` only when the browser should call an API on another origin. During local development, Vite proxies `/api` and `/health` to `http://127.0.0.1:8000`.

Docker is also supported:

```powershell
docker compose up --build
```

Then open `http://localhost:4173`.

## Reproduce the scientific data pipeline

The repository does not redistribute the upstream datasets. Download a named release from the [COCONUT download portal](https://coconut.naturalproducts.net/download) and [NPASS download page](https://bidd.group/NPASS/downloadnpass.html), retain the original file, and record that release with `--version`. The importer supports CSV, TSV, SDF, gzip, and ZIP streams. It does not load a whole release into memory.

Install the backend first (Python 3.11 or 3.12 is recommended):

```powershell
cd backend
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Fetch every ChEMBL page, normalize activities, and keep rejected rows separately:

```powershell
python -m app.pipeline fetch-chembl `
  --target CHEMBL335 --standard-type IC50 --page-size 1000 `
  --output ../data/raw/chembl335_ic50.raw.jsonl

python -m app.pipeline normalize-chembl `
  --input ../data/raw/chembl335_ic50.raw.jsonl `
  --output ../data/processed/chembl335_ic50.normalized.jsonl `
  --rejects ../data/processed/chembl335_ic50.rejects.jsonl
```

`fetch-chembl` follows ChEMBL's `page_meta.next` link until it is null. `--max-records` is available only as an explicit smoke-test cap; omission means no client-side 1,000-record ceiling. Retriable HTTP statuses (408, 425, 429, and 5xx), timeouts, and transport failures use bounded exponential retry with `Retry-After` support. Output is written through a `.part` file and accompanied by a timestamped SHA-256 manifest.

Build one cross-source natural-product index. These examples match the official NPASS 3.0 structure export and the tested COCONUT 2024-08 Drug Discovery ZIP; replace versions and member names when using another release:

```powershell
python -m app.pipeline index-natural-products `
  --source NPASS --version NPASS-3.0 `
  --input ../data/raw/NPASS3.0_naturalproducts_structure.txt `
  --output ../data/processed/natural_products.sqlite --batch-size 2000

python -m app.pipeline index-natural-products `
  --source COCONUT --version COCONUT-2024-08-DD `
  --input ../data/raw/COCONUT_2024_08_DrugDiscovery.tsv.zip `
  --member COCONUT_2024_08_DrugDiscovery.tsv `
  --output ../data/processed/natural_products.sqlite --batch-size 2000
```

For a quick real-file parser check, add `--max-records 1000`; omit it for a complete release. Multi-file ZIPs require an exact `--member`, preventing accidental ingestion of metadata files.

Verify normalized ChEMBL molecules against both natural-product sources:

```powershell
python -m app.pipeline crossmatch-natural-products `
  --input ../data/processed/chembl335_ic50.normalized.jsonl `
  --database ../data/processed/natural_products.sqlite `
  --output ../data/processed/chembl335_ic50.natural-products.jsonl `
  --rejects ../data/processed/chembl335_ic50.crossmatch-rejects.jsonl
```

### Normalization and provenance rules

- RDKit sanitization, cleanup, largest-fragment selection, and uncharging produce canonical isomeric SMILES, standard InChI, and standard InChIKey. Original SMILES and whether fragments were removed remain recorded.
- Invalid or missing structures are rejected with their source row; no structure or activity is imputed.
- Duplicate structures share one compound row, while each unique `(source, dataset_version, source_id)` keeps its own provenance record and raw source row.
- ChEMBL values are converted from M, mM, µM/uM, nM, or pM to molar units and pActivity only when the value is positive and the unit is supported. Relations such as `<` and `>` are retained as censored observations.
- Measurements receive an assay-context key derived from target, target organism, endpoint type, assay type, and BAO format. Downstream code must aggregate only within the same key; this pipeline never pools incompatible contexts.
- Reported pChEMBL values are retained and consistency-checked against the normalized value rather than silently overwritten.

Every JSONL and SQLite result remains traceable to its upstream record. The SQLite schema keeps normalized compounds, source-specific provenance, and rejected rows in separate tables.

### Verified integration snapshot (2026-10-20)

A live complete ChEMBL fetch for `CHEMBL335` / `IC50` traversed 38 API pages and returned 3,796 records (SHA-256 `ca16d730191ff6ea6b3adb07161f6891599ef3e2315b5cc21abc7a67b229536f`). Normalization accepted all 3,796 records into 8 distinct assay contexts, found 26 records without structures, and produced SHA-256 `442c53a321dea17c74533970a688efec5cddfd118b4c24ea4daa92e48d61288a`. These counts are a dated reproducibility observation, not fixed properties of the changing ChEMBL service.

The complete real NPASS 3.0 structure export was ingested: 203,390 rows read, 201,971 provenance records accepted, 198,028 unique standardized structures, 3,944 duplicate structures, 1 duplicate source record, and 1,418 invalid structures retained in the rejects table. The complete COCONUT 2024-08 Drug Discovery export added 695,128 provenance records from 695,133 rows; five invalid structures were rejected and 191,120 duplicate structures were linked to existing standardized compounds. The combined index contains 702,036 unique structures and preserves 897,099 source-specific provenance records. Counts depend on the selected upstream release. Automated tests use small deterministic fixtures and never substitute fabricated measurements for external data.

Train and scaffold-validate the molecular bioactivity models from normalized ChEMBL and structure-verified natural-product data:

```powershell
python -m app.pipeline train-activity `
  --activities ../data/processed/chembl335_ic50.normalized.jsonl `
  --natural-products ../data/processed/chembl335_ic50.natural-products.jsonl `
  --output-dir ../experiments/ml_activity_chembl335 `
  --folds 5 --bootstrap-iterations 1000 --seed 2026
```

The pipeline aggregates real compatible measurements by standardized InChIKey before splitting, generates Morgan radius-2 fingerprints plus ten RDKit descriptors, and compares a mean baseline, Random Forest, and XGBoost on identical Bemis–Murcko scaffold folds. A nested scaffold split calibrates 90% residual intervals. Tanimoto distance to training chemistry provides an explicit applicability-domain flag.

### Executed ML result (2026-10-20)

The real dataset contained 2,915 measurements aggregated into 2,444 structures across 1,086 scaffolds; 362 structures were NPASS-verified natural products. Five-fold scaffold-held-out results were:

| Model | MAE | RMSE | R² | NPASS MAE | 90% interval coverage |
|---|---:|---:|---:|---:|---:|
| Mean baseline | 0.738 | 0.925 | -0.006 | 0.662 | 0.900 |
| Random Forest | 0.455 | 0.597 | 0.581 | 0.396 | 0.905 |
| XGBoost | 0.432 | 0.575 | 0.612 | 0.368 | 0.924 |

XGBoost was selected by the prespecified outer-fold MAE. Its MAE 95% scaffold-cluster interval was 0.408–0.461. XGBoost-minus-Random-Forest MAE was -0.023 (paired scaffold-cluster 95% CI -0.038 to -0.010). The applicability diagnostic identified 315/2,444 outer-fold predictions as outside domain, where MAE increased from 0.420 to 0.511.

The uncertainty analysis has an important negative result: interval width had essentially no rank relationship with absolute error (Spearman -0.003). The intervals achieved aggregate coverage but should not be interpreted as molecule-specific calibrated confidence. See the [ML experiment protocol](docs/ml_experiment.md) and generated model card.

Run the leakage-resistant evidence-reversal experiment directly from normalized ChEMBL and natural-product crossmatch files:

```powershell
python -m app.pipeline evidence-reversal `
  --activities ../data/processed/chembl335_ic50.normalized.jsonl `
  --natural-products ../data/processed/chembl335_ic50.natural-products.jsonl `
  --output-dir ../experiments/evidence_reversal_chembl335 `
  --global-cutoff-year 2015 --k 10 --active-threshold 6.0 `
  --bootstrap-iterations 2000 --permutation-iterations 2000 --seed 2026
```

The command creates cohort CSVs, JSON metrics, input checksums, an SVG comparison, and a generated scientific report. Rankings use only earlier evidence; the latest publication year is held out independently for each molecule, with a global-2015 split as a stricter sensitivity analysis.

### Executed evidence-reversal result (2026-10-20)

The compatible context contained 2,915 measurements from 2,451 molecules. The latest-year benchmark contained 162 molecules, including 31 verified through NPASS. PhytoTruth increased all-compound nDCG@10 from 0.891 to 0.977, but the paired difference of 0.086 had a 95% molecule-bootstrap interval of 0.000 to 0.265 and was not significant by paired randomization (`p=0.3458`). In the NPASS subgroup, nDCG@10 decreased from 0.891 to 0.765.

The ablation is also negative for the reliability additions: median potency alone matched the full heuristic's top-10 metrics and produced better global ordering. Full-minus-median Kendall tau-b was -0.035 (95% CI -0.079 to -0.003; `p=0.0400`). Therefore this experiment does **not** validate an evidence-aware advantage; its all-compound point-estimate gain is attributable to robust central potency rather than demonstrated benefit from the reliability components. See the [experiment protocol](docs/experiment.md) for the complete assumptions and limitations.

### Robust target-held-out evaluation (2026-10-26)

The follow-up runner discovered 30 local target datasets and retained 20 using a coverage-only rule. Before fitting, 12 complete targets were assigned to development and 8 to locked evaluation. Within every target, each molecule's latest publication year was held out. Ridge strength and reliability feature family were selected exclusively by leave-one-development-target-out validation.

The selected `reliability:ridge:10` candidate did not pass the development gate. On 1,235 held-out molecules across the 8 evaluation targets, macro nDCG@10 changed from 0.9341 to 0.9348: difference +0.0007, 95% target-bootstrap CI -0.0053 to +0.0084, exact target sign-flip `p=0.890625`, Holm-adjusted `p=1.0`, paired Hedges g=0.060, with 2/0/6 target wins/ties/losses. Kendall tau changed by -0.0088 (95% CI -0.0160 to -0.0014; exact `p=0.0625`; Holm-adjusted `p=0.25`). Median potency therefore remains the deployed guardrail.

Eight independent evaluation targets permit a minimum two-sided exact p-value of only 0.0078125, so this run could not honestly reach p<=0.001. The generated artifacts preserve that power limitation and all negative ablations. See the [robust multi-target protocol](docs/robust-multitarget-evaluation.md) and `experiments/robust_multitarget_evaluation/report.md`.

### Study-aware Evidence Intelligence

The decision-support compiler uses the executed development-only probabilistic benchmark and selects the method minimizing worst NLL rank across study-grouped and forward-temporal validation. Student-t regression ranks first and second respectively. In the 25-event forward-time evaluation it achieved MAE 0.517 and 84% coverage for nominal 90% intervals. The derived ten-fold-disagreement probability had AUROC 0.826, but only three disagreements were available and the study-cluster 95% CI was 0.375–0.958; it is therefore not a validated classifier.

Uncertainty-based abstention was more promising: retaining the 60% narrowest intervals reduced retrospective MAE to 0.317 and raised observed 90% coverage to 100%. The Study Planner turns this distribution into outcome-blind experiment recommendations using a model-based expected-information-gain score. These recommendations prioritize independent information collection and are not efficacy or success predictions. See [the Study-aware Evidence Intelligence protocol](docs/study-aware-evidence-intelligence.md).

### Evidence Reliability investigation

The follow-up investigation uses the same strict latest-publication holdout, then creates a deterministic molecule-level development/evaluation partition using only molecule IDs and NPASS membership. A prespecified development-only gate determines whether the existing evidence adjustment may replace median potency; evaluation outcomes never enter that decision. It compares maximum, mean, most-recent, lower-confidence-bound, and median potency baselines, quantifies component predictive value, and writes compound-level harmful reversal attribution.

```powershell
cd backend
python -m app.reliability_investigation `
  --activities "..\..\work\chembl335_ic50.normalized.jsonl" `
  --natural-products "..\..\work\chembl335_ic50.npass.jsonl" `
  --output-dir "..\experiments\evidence_reliability_chembl335" `
  --development-fraction 0.6 --k 10 --active-threshold 6.0 `
  --iterations 2000 --seed 2026
```

The resulting `report.md`, `results.json`, metric/component CSVs, comparison SVG, compound error tables, input hashes, and exact command form a reproducible audit trail. Reliability values remain heuristic descriptors, not calibrated probabilities. The current natural-product subgroup is verified through NPASS only; COCONUT did not contribute records to this target artifact.

### Multi-target scientific benchmark

The expanded benchmark prespecifies five human single-protein targets before examining ranking outcomes: acetylcholinesterase (`CHEMBL220`), EGFR (`CHEMBL203`), COX-2 (`CHEMBL230`), DPP-4 (`CHEMBL284`), and PTP1B (`CHEMBL335`). Each target must have at least 30 molecules with strict latest-publication holdout evidence and at least five structure-verified natural products in its largest compatible assay context. Context selection maximizes molecule coverage and never uses activity values.

Fetch and normalize each target with the same uncapped pagination and normalization commands shown above, then crossmatch every normalized file against the combined NPASS/COCONUT SQLite index. Record the resulting paths and authoritative ChEMBL target names in `experiments/multitarget_targets.json`, then run:

```powershell
cd backend
python -m app.pipeline multi-target-benchmark `
  --specs ../experiments/multitarget_targets.json `
  --natural-product-database ../data/processed/natural_products.sqlite `
  --output-dir ../experiments/multi_target_benchmark `
  --development-fraction 0.6 --k 10 --active-threshold 6.0 `
  --iterations 1000 --seed 2026
```

The benchmark performs publication-time separation first and a molecule-level development/evaluation partition second. Across-target intervals use hierarchical target-then-molecule resampling. Median potency remains the ranking guardrail unless an evidence-aware alternative passes the prespecified development gate. See [the multi-target protocol](docs/multitarget-benchmark.md); generated artifacts record input hashes, exclusions, target-specific metrics, natural-product source coverage, contradictions, transferability comparisons, and compound-level harmful reversals.

### Evidence transferability and assay disagreement

The follow-up treats later-study disagreement as its own outcome instead of retuning ranking weights. A compound event is eligible only when its latest publication-year exact-IC50 evidence comes from document- and assay-disjoint records in the same compatible target context as its earlier evidence. The prespecified primary outcome is a ten-fold or larger absolute potency shift.

```powershell
cd backend
python -m app.pipeline assay-disagreement `
  --specs ../experiments/multitarget_targets.json `
  --output-dir ../experiments/assay_disagreement `
  --bootstrap-iterations 1000 `
  --prediction-interval-iterations 200 `
  --seed 2026
```

Development outcomes end in 2019 and untouched evaluation begins in 2020; molecules are disjoint across partitions. Predictors use only RDKit descriptors and evidence known before the held-out study. NPASS/COCONUT membership is reserved for subgroup reporting. Inference resamples publication clusters within target because many compounds can share the same experimental study. See [the assay-disagreement protocol](docs/assay-disagreement.md).

The executed experiment retained 531 independent later-study events after quality filtering: 212 development events (50 disagreements) and 319 untouched evaluation events (169 disagreements). The primary full model had AUROC 0.642 (study-clustered 95% CI 0.414–0.686), PR-AUC 0.609 (95% CI 0.303–0.654), and Brier score 0.296. Its PR-AUC gain over the constant development-prevalence baseline was 0.079 (paired study-clustered 95% CI -0.042 to 0.120). It therefore failed the prespecified support rule and is **not a validated disagreement predictor**. The natural-product evaluation was exploratory (60 compounds, 19 disagreements; AUROC 0.483, PR-AUC 0.336). One EGFR publication contained 131 compounds and 102 disagreements, demonstrating substantial study-level clustering and protocol-shift risk.

### Study-grouped generalization investigation

The robustness follow-up uses leave-one-publication-out predictions for every development compound and cross-fits calibration by publication. Model selection minimizes equal-study-weighted Brier score using development only. This prevents a large compound series from determining the selected method.

```powershell
cd backend
python -m app.pipeline assay-disagreement-generalization `
  --specs ../experiments/multitarget_targets.json `
  --output-dir ../experiments/assay_disagreement_generalization `
  --bootstrap-iterations 1000 `
  --prediction-interval-iterations 200 `
  --seed 2026
```

The negative result is scientifically important. The constant development-prevalence fallback had the best study-balanced development Brier score (`0.196`); compact logistic (`0.201`), shallow gradient boosting (`0.202`), history-only logistic (`0.204`), calibrated full logistic (`0.208`), and the existing uncalibrated full logistic (`0.227`) were worse. The fallback was therefore frozen before locked evaluation. It has AUROC `0.500`, PR-AUC `0.530`, and Brier `0.336`; it is less discriminative than the prior model's point estimates but avoids promoting a signal that failed study-grouped development validation. The status remains **not supported for reliable prediction**. See [the study-grouped protocol](docs/assay-disagreement-generalization.md).

### Assay protocol metadata enrichment

The protocol-enrichment phase retrieves the exact ChEMBL assay records referenced by the frozen 531-event dataset. It preserves each raw response, request URL, retrieval timestamp, record checksum, assay/document/target identifiers, controlled assay fields, parameters, variant annotations, and verbatim description. A versioned rule engine extracts only literal mentions of methodology, detection, biological system, endpoint, substrate, temperature, pH, and duration; an empty feature means *unreported or not detected*, not experimental absence.

```powershell
cd backend
python -m app.pipeline fetch-assay-metadata `
  --events ../experiments/assay_disagreement_generalization/dataset.csv `
  --output ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl `
  --batch-size 50 --max-retries 5

python -m app.pipeline enrich-assay-protocols `
  --events ../experiments/assay_disagreement_generalization/dataset.csv `
  --assays ../experiments/assay_protocol_enrichment/assay_metadata.raw.jsonl `
  --output-dir ../experiments/assay_protocol_enrichment `
  --bootstrap-iterations 1000 --seed 2026
```

All 1,050 requested assays were retrieved in the executed snapshot. Prior protocol fields are prediction-time eligible; held-out protocol fields and protocol-difference flags are explicitly `post_outcome_explanatory_only`. The analysis is descriptive and study-clustered and does not train or retrain a model. See [the assay protocol enrichment protocol](docs/assay-protocol-enrichment.md) and the generated scientific report under `experiments/assay_protocol_enrichment/`.

### Independent frozen-model validation

The independent replication uses five previously unseen human targets and publications from 2020 onward. The target list and sufficiency criteria were frozen before constructing disagreement outcomes. Every compound, activity, assay, and publication found anywhere in the previous five normalized target snapshots is blocked; compounds repeated across new targets are also removed. Both existing disagreement artifacts are checksum-verified and scored without fitting, calibration, or tuning.

```powershell
cd backend
python -m app.pipeline independent-assay-validation `
  --specs ../experiments/independent_validation_targets.json `
  --protocol ../experiments/independent_validation_protocol.json `
  --previous-specs ../experiments/multitarget_targets.json `
  --previous-events ../experiments/assay_disagreement_generalization/dataset.csv `
  --full-model ../experiments/assay_disagreement/assay_disagreement_model.joblib `
  --selected-model ../experiments/assay_disagreement_generalization/study_generalized_model.joblib `
  --output-dir ../experiments/independent_assay_validation `
  --bootstrap-iterations 2000 --seed 2026
```

The executed dataset retained 90 compounds from 37 independent held-out publications, but only seven disagreement outcomes and five verified natural products. It therefore failed the prespecified minimum size and cannot establish reliable independent performance. The frozen full model also failed the performance rule (AUROC `0.494`, PR-AUC `0.093`, Brier `0.421`). See [the independent-validation protocol](docs/independent-assay-validation.md) and the generated scientific report.

## Semantic knowledge graph and virtual screening

Build the evidence graph from normalized target activities and the structure-provenance index:

```powershell
cd backend
python -m app.pipeline build-knowledge-graph `
  --specs ../examples/knowledge_graph_targets.example.json `
  --natural-product-database D:\data\natural_products.sqlite `
  --output ../artifacts/knowledge_graph.sqlite
```

Run target- or disease-specific batch screening. The index is streamed, so the complete
NPASS/COCONUT collection can be evaluated without retaining it in memory:

```powershell
python -m app.pipeline screen-natural-products `
  --specs ../examples/screening_spec.example.json `
  --natural-product-database D:\data\natural_products.sqlite `
  --knowledge-graph ../artifacts/knowledge_graph.sqlite `
  --output-dir ../artifacts/ptp1b_screen `
  --top-k 250
```

The graph stores typed nodes and evidence-bearing edges for organisms, compounds, natural-
product records, targets, diseases, assays and publications. Screening writes ranked CSV and
JSONL candidates, a checksum-bearing manifest, full per-target predictions, residual intervals,
applicability-domain results, physicochemical descriptors, source provenance, and structural
alerts. These artifacts prioritize laboratory follow-up; they are not efficacy, safety, or
clinical claims. See [the graph and screening guide](docs/knowledge-graph-screening.md).

## Verify

```powershell
cd backend
pytest -q
```

The test suite also covers temporal splitting, activity/document leakage prevention, deterministic ranking, molecular feature generation, scaffold isolation, deterministic model fitting, native model loading and checksums, API input/error handling, deterministic predictions, applicability-domain warnings, held-out-evidence isolation, assay-metadata batching, literal protocol extraction, missingness preservation, provenance joins, feature-timing labels, clustered-analysis reproducibility, and duplicate-metadata rejection.

## Scientific sources

- [ChEMBL REST API documentation](https://www.ebi.ac.uk/chembl/api/data/docs)
- [COCONUT 2.0 download portal](https://coconut.naturalproducts.net/download) and [COCONUT 2.0 paper](https://doi.org/10.1093/nar/gkae1063)
- [NPASS official downloads](https://bidd.group/NPASS/downloadnpass.html) and [NPASS 2023 update](https://doi.org/10.1093/nar/gkac1069)
- Morgan, H. L. (1965), connectivity-based molecular fingerprints; RDKit implementation is used here.
- Bemis & Murcko (1996), framework-based scaffold grouping.

See [architecture](docs/architecture.md) and [experiment protocol](docs/experiment.md) for assumptions, leakage controls and evaluation design.

## Integrity checklist

- No generated assay values or fabricated citations.
- Every bundled measurement preserves a ChEMBL activity, assay and document ID.
- Measured evidence and model predictions are separate data types.
- Reliability is called a heuristic score, not probability or truth.
- Natural-product provenance is not inferred from a compound name.
- Evaluation is against held-out evidence, not the ranking heuristic itself.

## Current data limitations

- Upstream records with absent or invalid structures remain explicit rejects; the pipeline does not repair them heuristically.
- The compact COCONUT Drug Discovery export contains structures and derived scores but not all organism/reference fields available in richer exports. Use a fuller official export when those fields are required.
- The NPASS structure file verifies NPASS membership. NPASS organism and activity downloads are separate files and are not joined in this phase; their evidence must not be attributed to a compound until a future keyed join is implemented and tested.
- Assay-context keys prevent accidental pooling but do not by themselves prove two assays are biologically equivalent.

# Independent assay-disagreement dataset expansion

The `assay-disagreement-dataset-expansion-v1` phase adds a label-quarantined,
strictly independent cohort from 15 prespecified human single-protein targets.
Target selection and the 1.0-pChEMBL disagreement threshold were frozen before
outcomes were constructed. The pipeline uses exact IC50 values only, selects one
target-specific assay context without activity outcomes, requires earlier and
latest-year evidence to have disjoint ChEMBL assay and document identifiers, and
removes every compound, activity, assay, publication, and target exposed in the
two earlier assay-disagreement cohorts.

Natural-product status is an exact standardized-InChIKey match against the local
NPASS + COCONUT provenance index. It is database membership, not evidence of
efficacy or clinical safety. Missing protocol fields remain unknown.

From `backend`, reproduce the snapshots and expansion as follows (repeat the
three data commands for every target in
`experiments/assay_expansion_targets.json`):

```powershell
python -m app.pipeline fetch-chembl --target CHEMBL1951 --standard-type IC50 --output ..\data\assay_disagreement_expansion\CHEMBL1951.raw.jsonl --page-size 1000
python -m app.pipeline normalize-chembl --input ..\data\assay_disagreement_expansion\CHEMBL1951.raw.jsonl --output ..\data\assay_disagreement_expansion\CHEMBL1951.normalized.jsonl --rejects ..\data\assay_disagreement_expansion\CHEMBL1951.normalize-rejects.jsonl
python -m app.pipeline crossmatch-natural-products --input ..\data\assay_disagreement_expansion\CHEMBL1951.normalized.jsonl --database PATH_TO_NPASS_COCONUT.sqlite --output ..\data\assay_disagreement_expansion\CHEMBL1951.natural.jsonl --rejects ..\data\assay_disagreement_expansion\CHEMBL1951.natural-rejects.jsonl

python -m app.dataset_expansion prepare --specs ..\experiments\assay_expansion_targets.json --data-root ..\data\assay_disagreement_expansion --previous-specs ..\experiments\multitarget_targets.json --previous-events ..\experiments\assay_disagreement_generalization\dataset.csv --previous-specs ..\experiments\independent_validation_targets.json --previous-events ..\experiments\independent_assay_validation\independent_dataset.csv --output-dir ..\experiments\assay_disagreement_expansion
python -m app.assay_protocol fetch --events ..\experiments\assay_disagreement_expansion\expansion_events.csv --output ..\experiments\assay_disagreement_expansion\assays.raw.jsonl --batch-size 50
python -m app.dataset_expansion finalize --events ..\experiments\assay_disagreement_expansion\expansion_events.csv --measurements ..\experiments\assay_disagreement_expansion\used_measurements.csv --assays ..\experiments\assay_disagreement_expansion\assays.raw.jsonl --output-dir ..\experiments\assay_disagreement_expansion\final
pytest -q
```

The source snapshots have checksum manifests. The experiment folder contains the
event table, all 295 used measurements as both a compact CSV and complete
normalized JSONL records, exact NPASS/COCONUT provenance, explicit exclusions,
196 raw assay records, extracted protocol fields, positive-case studies, coverage
statistics, a scientific report, and output checksums. No model is loaded,
trained, tuned, or evaluated by this phase.

## Expanded frozen-model validation

The 94-event expansion was subsequently used as a locked test set for the two
existing assay-disagreement artifacts. The audit found no overlap with either
earlier ledger at compound, target, activity, assay, or publication level, so all
94 events remained eligible. Model files were checksum-verified and were not
fitted, recalibrated, tuned, or selected.

The frozen full model did not generalize: overall AUROC was `0.390` (study-
clustered 95% CI `0.240–0.540`), PR-AUC `0.225`, Precision@10 `0.300`, and Brier
`0.244`. The frozen prevalence fallback had Brier `0.191`. Among 20 verified
natural products, frozen-full AUROC was `0.209` (`0.107–0.338`), PR-AUC `0.261`,
Precision@10 `0.100`, and Brier `0.346`; all seven disagreements were below the
fixed 0.5 descriptive threshold. This is a negative result, and no winning model
is selected. See [the complete expanded validation report](docs/expanded-assay-validation.md).

## Development-only disagreement-methodology redesign

The redesign intentionally does **not** load either independent validation cohort. It uses only the original pre-2020 development partition, changes the primary outcome to continuous absolute pChEMBL shift, and retains the 1.0-pChEMBL binary label as a secondary view. Validation is molecule- and publication-grouped, forward-temporal, and leave-one-target-out. Candidate hyperparameters and readiness rules are recorded in `experiments/methodology_redesign_protocol.json`.

```powershell
cd backend
python -m app.methodology_redesign `
  --dataset ../experiments/assay_protocol_enrichment/enriched_assay_disagreement_events.csv `
  --output-dir ../experiments/assay_disagreement_methodology_redesign `
  --bootstrap-iterations 1000 --seed 2026
python -m pytest tests/test_methodology_redesign.py -q
python -m pytest -q
```

Outputs include fold-level predictions, clustered confidence intervals, paired protocol-feature ablations, natural-product error analysis, an SVG comparison, candidate artifacts, and a generated scientific report. Candidate models are development artifacts and are explicitly marked as not independently validated; they do not replace the frozen production or prior research models.

## PhytoTRUST-M0 study-aware feasibility experiment

PhytoTRUST-M0 is the stage-gated hierarchical baseline specified in
`docs/phytotrust-research-methodology.md`. It predicts a future study's pChEMBL
measurement as the historical compound-target median plus a target correction,
a partially pooled target-specific study random intercept, and Student-t residual
noise. For blind transfer, the unseen study effect is integrated out; no outcome
from the held-out study is used to estimate it.

Run the development-only experiment from `backend/`:

```powershell
python -m app.phytotrust_m0 `
  --dataset ..\experiments\assay_protocol_enrichment\enriched_assay_disagreement_events.csv `
  --output-dir ..\experiments\phytotrust_m0 `
  --bootstrap-iterations 1000 `
  --seed 2026
```

Run the M0 tests or the complete backend suite:

```powershell
python -m pytest tests\test_phytotrust_m0.py -q
python -m pytest -q
```

The experiment excludes the locked 2020+ evaluation partition and both previously
examined independent cohorts. Its outputs are development-only feasibility
evidence, not independent validation or proof of methodological novelty. See
`experiments/phytotrust_m0/SCIENTIFIC_REPORT.md` for exact results,
study-cluster bootstrap comparisons, identifiability diagnostics, and limitations.

## PhytoTRUST study-connected dataset

The `phytotrust-connected-dataset-v1` protocol freezes ten previously unused human
single-protein targets and identifier-only connectivity criteria before target data
are evaluated. Target/context selection never uses IC50 values. Records sharing a
compound, activity, assay, or publication with any previously locked evaluation
cohort are excluded.

From `backend/`, retrieve each target listed in
`experiments/phytotrust_connected_dataset_protocol.json`, then normalize and
crossmatch it against the existing NPASS/COCONUT index:

```powershell
python -m app.pipeline fetch-chembl --target CHEMBL222 --standard-type IC50 `
  --output ..\data\phytotrust_connected\CHEMBL222.raw.jsonl --page-size 1000
python -m app.pipeline normalize-chembl `
  --input ..\data\phytotrust_connected\CHEMBL222.raw.jsonl `
  --output ..\data\phytotrust_connected\CHEMBL222.normalized.jsonl `
  --rejects ..\data\phytotrust_connected\CHEMBL222.normalize-rejects.jsonl
python -m app.pipeline crossmatch-natural-products `
  --input ..\data\phytotrust_connected\CHEMBL222.normalized.jsonl `
  --database PATH_TO_NPASS_COCONUT.sqlite `
  --output ..\data\phytotrust_connected\CHEMBL222.natural.jsonl `
  --rejects ..\data\phytotrust_connected\CHEMBL222.natural-rejects.jsonl
```

After all frozen candidates are present:

```powershell
python -m app.phytotrust_connected_dataset fetch-target-metadata `
  --protocol ..\experiments\phytotrust_connected_dataset_protocol.json `
  --output ..\data\phytotrust_connected\target_metadata.json
python -m app.phytotrust_connected_dataset build `
  --protocol ..\experiments\phytotrust_connected_dataset_protocol.json `
  --target-metadata ..\data\phytotrust_connected\target_metadata.json `
  --data-root ..\data\phytotrust_connected `
  --experiments-root ..\experiments `
  --output-dir ..\experiments\phytotrust_connected_dataset
python -m app.phytotrust_connected_dataset validate `
  --output-dir ..\experiments\phytotrust_connected_dataset `
  --experiments-root ..\experiments
python -m pytest tests\test_phytotrust_connected_dataset.py -q
```

Development and future evaluation are split by complete compound-publication graph
component, so compounds, publications, assays, and activities cannot cross the
partition boundary. Pre-2020 evidence belonging to future components is stored as
evaluation history, while 2020+ values are written only to the explicitly
quarantined label file. Passing the structural graph gate is necessary but does not
prove that a hierarchical variance component is statistically identifiable.

## PhytoTRUST-M0 on the connected CHEMBL222 development graph

This re-evaluation consumes only
`experiments/phytotrust_connected_dataset/development_measurements.csv`. The loader
fails closed for future-evaluation or quarantined artifact paths. Compatible records
are aggregated to compound-publication medians, while every activity and assay ID is
retained for fold-level leakage checks.

Run the frozen development-only experiment from `backend/`:

```powershell
python -m app.phytotrust_m0_connected `
  --dataset ..\experiments\phytotrust_connected_dataset\development_measurements.csv `
  --output-dir ..\experiments\phytotrust_m0_connected `
  --bootstrap-iterations 1000 `
  --parameter-bootstrap-iterations 200 `
  --seed 2026
python -m pytest tests\test_phytotrust_m0_connected.py -q
python -m pytest -q
```

The study-grouped analysis is a retrospective publication-domain holdout; only the
forward-temporal analysis is prospective with respect to publication year. Earlier
or otherwise allowed measurements of the same compound are part of the declared
transfer estimand, but the held-out publication, activity, assay, and outcome are
excluded. Results are development-only feasibility evidence and must not be treated
as independent validation, clinical evidence, or proof of methodological novelty.

## PhytoTRUST-M1 evidence-conditioned transferability

M1 extends the connected M0 feasibility model with a regularized historical-
evidence correction and an evidence-conditioned Student-t scale. Only protocol
metadata attached to allowed historical measurements are predictors; held-out
publication protocol fields remain unavailable. The frozen protocol is
`experiments/phytotrust_m1_protocol.json`, and the mathematical specification is
documented in `docs/phytotrust-m1-methodology.md`.

Retrieve and extract the 67 development-only ChEMBL assay records from `backend/`:

```powershell
python -m app.phytotrust_m1 fetch-assays `
  --dataset ..\experiments\phytotrust_connected_dataset\development_measurements.csv `
  --raw-output ..\data\phytotrust_m1\CHEMBL222.assays.raw.jsonl `
  --feature-output ..\data\phytotrust_m1\CHEMBL222.assay_protocol_features.csv `
  --batch-size 50 --max-retries 5
```

Run the frozen development-only comparison:

```powershell
python -m app.phytotrust_m1 run `
  --dataset ..\experiments\phytotrust_connected_dataset\development_measurements.csv `
  --protocol-features ..\data\phytotrust_m1\CHEMBL222.assay_protocol_features.csv `
  --output-dir ..\experiments\phytotrust_m1 `
  --bootstrap-iterations 1000 --seed 2026
python -m pytest tests\test_phytotrust_m1.py -q
python -m pytest -q
```

Outputs include paired fold predictions, baseline and ablation tables,
publication-cluster confidence intervals, uncertainty components, two SVG figures,
model configuration/parameters, source hashes, and a generated scientific report.
The support rule is evaluated before any quarantined final-test access. A failed rule
is a stopping result, not permission to tune on the final cohort.

