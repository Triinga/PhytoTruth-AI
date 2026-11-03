# Semantic knowledge graph and natural-product screening

PhytoTruth can build a portable evidence graph and run a disk-backed, bounded-memory
multi-target screen over the existing NPASS/COCONUT index. These are research
prioritization tools. They do not establish efficacy, safety, causality, or clinical utility.

## Graph semantics

The SQLite artifact contains typed nodes for `organism`, `compound`,
`natural_product_record`, `target`, `disease`, `assay`, and `publication`. Edges include:

- `REPORTED_SOURCE_OF`: an upstream natural-product record reports an organism as a source;
- `VERIFIED_BY_NATURAL_PRODUCT_SOURCE`: a standardized InChIKey matches NPASS/COCONUT;
- `MEASURED_AGAINST`, `TESTED_IN`, `MEASURES`, and `REPORTS`: ChEMBL experimental evidence;
- `ASSOCIATED_WITH_DISEASE`: an explicitly curated relationship from the input specification;
- `PREDICTED_AGAINST` and `PRIORITIZED_FOR_DISEASE`: outputs added by a screening run.

Every edge stores `evidence_source`, `evidence_id`, and relationship properties. The graph
does not infer that an organism is a plant; upstream organism strings are retained with
`plant_status=not_taxonomically_asserted` until a taxonomy source is joined.

Example graph specification:

```json
{
  "targets": [
    {
      "target_chembl_id": "CHEMBL335",
      "target_name": "Protein-tyrosine phosphatase 1B",
      "target_organism": "Homo sapiens",
      "activities": "../data/processed/chembl335_ic50.normalized.jsonl",
      "diseases": [
        {"id": "MONDO:0005148", "name": "type 2 diabetes mellitus"}
      ]
    }
  ]
}
```

Build it:

```powershell
cd backend
python -m app.pipeline build-knowledge-graph `
  --specs ../examples/knowledge_graph_targets.example.json `
  --natural-product-database D:\data\natural_products.sqlite `
  --output ../artifacts/knowledge_graph.sqlite
```

The command also writes `knowledge_graph.sqlite.manifest.json` with checksums, counts,
versions, paths, and semantics.

## Multi-target screening

Each target must point to a checksum-validated native inference model produced by
`train-activity`. The screening score is a transparent weighted aggregation of predicted
pChEMBL, applicability-domain membership, nearest-training Tanimoto similarity, interval
width, and structural-alert burden. Physicochemical filters are configurable. PAINS/Brenk
matches are reported as research alerts and are not represented as toxicity predictions.

```powershell
python -m app.pipeline screen-natural-products `
  --specs ../examples/screening_spec.example.json `
  --natural-product-database D:\data\natural_products.sqlite `
  --knowledge-graph ../artifacts/knowledge_graph.sqlite `
  --output-dir ../artifacts/ptp1b_screen `
  --top-k 250
```

Use `--source NPASS` or `--source COCONUT` to restrict provenance. Use
`--max-compounds 1000` only for smoke tests. Without it, the index is streamed completely;
only the top K candidates remain in memory.

Outputs:

- `ranked_candidates.jsonl`: complete nested predictions and provenance for the retained set;
- `ranked_candidates.csv`: analyst-friendly export;
- `screening_manifest.json`: hashes, filters, counts, model targets, scoring semantics, and limitations;
- optional graph edges connecting retained compounds to targets and the disease.

## API exposure

The API serves precomputed artifacts and never starts an expensive screen inside a request:

```powershell
$env:PHYTOTRUTH_KNOWLEDGE_GRAPH = "..\artifacts\knowledge_graph.sqlite"
$env:PHYTOTRUTH_SCREENING_RESULTS = "..\artifacts\ptp1b_screen\ranked_candidates.jsonl"
uvicorn app.main:app --reload
```

- `GET /api/knowledge-graph/{node_type}/{external_id}` returns a bounded evidence neighborhood.
- `GET /api/screening/candidates?limit=25` returns precomputed ranked candidates.

