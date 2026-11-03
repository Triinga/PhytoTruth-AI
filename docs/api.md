# Bioactivity prediction API

`POST /api/predict` runs the saved CHEMBL335 XGBoost model with the exact RDKit feature configuration used during training: a 2,048-bit radius-2 Morgan fingerprint followed by ten ordered molecular descriptors. The endpoint is for research prioritization only. Its output is a target-specific pChEMBL estimate, not confirmed activity, drug efficacy, clinical benefit, or safety.

## Model artifacts

Runtime loading uses XGBoost's supported native JSON format rather than unpickling an XGBoost estimator. The applicability-domain reference fingerprints are stored in a NumPy archive loaded with `allow_pickle=False`. SHA-256 checksums, feature dimensions, and artifact versions are verified before inference.

The checked-in inference set is:

```text
experiments/ml_activity_chembl335/
  inference_config.json
  inference_manifest.json
  model.xgb.json
  reference_fingerprints.npz
```

To reproducibly convert the existing `best_model.joblib` without retraining:

```powershell
cd backend
python -m app.model_export `
  --source ../experiments/ml_activity_chembl335/best_model.joblib `
  --output-dir ../experiments/ml_activity_chembl335
```

The legacy bundle is read only by this conversion command. API requests load the native artifacts. Keep the conversion environment version-consistent with the training environment; the exported XGBoost version and source hash are recorded in `inference_config.json`.

Set `PHYTOTRUTH_MODEL_DIR` to use another artifact directory:

```powershell
$env:PHYTOTRUTH_MODEL_DIR = "C:\models\phytotruth"
uvicorn app.main:app --host 127.0.0.1 --port 8000
```

The first request validates and loads the artifact. The process caches the resulting predictor for subsequent requests.

## Request

```bash
curl -X POST http://127.0.0.1:8000/api/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'
```

PowerShell:

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/api/predict `
  -ContentType application/json `
  -Body '{"smiles":"CCO"}'
```

## Example response

This response was produced by model version `chembl335-xgboost-c20d771793ac`; the values are model outputs, not measured ethanol bioactivity:

```json
{
  "canonical_smiles": "CCO",
  "predicted_pchembl": 4.849597454071045,
  "prediction_interval": {
    "lower": 2.988195328965653,
    "upper": 6.710999579176437,
    "confidence_level": 0.9,
    "method": "cross-validated absolute-residual conformal heuristic"
  },
  "applicability_domain": {
    "inside": false,
    "status": "OUTSIDE_DOMAIN_DO_NOT_RELY",
    "nearest_training_tanimoto": 0.21052631578947367,
    "threshold": 0.4225352112676056,
    "warning": "Outside applicability domain; interval calibration is not guaranteed."
  },
  "descriptors": {
    "molecular_weight": 46.069,
    "logp": -0.001,
    "tpsa": 20.23,
    "hbd": 1.0,
    "hba": 1.0,
    "rotatable_bonds": 0.0,
    "fraction_csp3": 1.0,
    "ring_count": 0.0,
    "heavy_atoms": 3.0,
    "formal_charge": 0.0
  },
  "model": {
    "name": "xgboost",
    "version": "chembl335-xgboost-c20d771793ac",
    "format": "xgboost-json",
    "target": {
      "assay_type": "B",
      "bao_format": "BAO_0000357",
      "standard_type": "IC50",
      "target_chembl_id": "CHEMBL335",
      "target_organism": "Homo sapiens"
    },
    "training_structures": 2444,
    "training_scaffolds": 1086
  },
  "kind": "predicted",
  "limitations": [
    "OUTSIDE APPLICABILITY DOMAIN: do not treat this prediction as reliable.",
    "Research prioritization only; not confirmed experimental activity.",
    "Not evidence of drug efficacy, clinical benefit, or safety.",
    "Intervals are retrospective cross-validated residual diagnostics, not guaranteed prospective coverage.",
    "Predictions outside the applicability domain must not be treated as reliable."
  ]
}
```

## Errors

- `422 Unprocessable Entity`: blank, oversized, or RDKit-invalid SMILES. Invalid structures receive the machine-readable code `invalid_smiles` when request parsing succeeds.
- `503 Service Unavailable`: the model directory, configuration, native model, or fingerprint archive is missing, corrupt, checksum-invalid, or incompatible.
- `500 Internal Server Error`: an unexpected failure occurred after a valid artifact and molecule were loaded.

## Interpretation and limitations

- The 90% interval uses the stored 90th percentile of outer scaffold-held-out absolute residuals. It is a retrospective diagnostic, not a probability that this individual prediction is correct.
- For structures below the training-set Tanimoto threshold, the interval is widened heuristically and the result is explicitly marked `OUTSIDE_DOMAIN_DO_NOT_RELY`. Calibration is not guaranteed there.
- The training experiment found essentially no rank relationship between interval width and absolute error. Do not interpret interval width as molecule-specific calibrated confidence.
- The model estimates the compatible CHEMBL335 IC50 context recorded in the response. It must not be generalized to other targets, assay contexts, organisms, efficacy, dosing, toxicity, or clinical safety.

Interactive OpenAPI documentation is available at `http://127.0.0.1:8000/docs` while the API is running.
