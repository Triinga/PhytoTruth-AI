# Multi-target evidence benchmark protocol

## Scientific question

Does assay reliability information improve prioritization of compounds whose later, independently published IC50 evidence is unknown? The benchmark tests this across targets; it does not assess clinical efficacy or safety.

## Prespecified targets

- `CHEMBL220`: acetylcholinesterase
- `CHEMBL203`: epidermal growth factor receptor
- `CHEMBL230`: prostaglandin G/H synthase 2 (COX-2)
- `CHEMBL284`: dipeptidyl peptidase 4
- `CHEMBL335`: tyrosine-protein phosphatase non-receptor type 1 (PTP1B)

All are human single-protein ChEMBL targets. They were chosen for biological diversity and adequate activity coverage before benchmark outcomes were examined.

## Data integrity

ChEMBL IC50 pages are followed until `page_meta.next` is null. Exact positive concentration measurements in supported units are converted to molar concentration and pActivity. Activity IDs are unique, and raw values, units, reported pChEMBL, validity annotations, duplicate flags, document IDs, and assay metadata are retained. Censored, malformed, structurally invalid, or incompatible records are never imputed or silently combined.

NPASS 3.0 and COCONUT 2024-08 structures are standardized with RDKit cleanup, largest-fragment selection, uncharging, canonical isomeric SMILES, standard InChI, and InChIKey. A compound is a verified natural product only through standardized structure identity. Source-specific IDs and raw rows remain in the provenance table when sources overlap.

## Leakage controls

1. Within each target, the compatible assay context is chosen by distinct-molecule coverage without using activity values.
2. For each molecule, all evidence from its latest publication year is held out. Training publications must be strictly earlier and document-disjoint.
3. Molecules are deterministically hashed into 60% development and 40% evaluation partitions, stratified only by natural-product membership.
4. The current PhytoTruth heuristic and nonparametric 95% lower-confidence-bound alternative can replace median potency only if development macro nDCG@10 and Kendall tau-b both have hierarchical-bootstrap interval lower bounds above zero and paired score-swap p-values below 0.05.
5. The method is frozen before evaluation outcomes are summarized.

The evaluation is untouched by this revision procedure but is not a prospective external dataset.

## Evaluation

Per target and natural-product subgroup, the benchmark reports nDCG@10, precision@10, enrichment factor, Kendall tau-b, every existing ablation, and strong potency baselines. Across-target comparisons use equal-weight macro averages. Confidence intervals resample targets first and molecules within targets second; randomization swaps paired method scores within targets.

Compound-level error tables identify harmful and beneficial pairwise ranking reversals and decompose the exact contribution from assay support, inconsistency, composite reliability, and bootstrap uncertainty.

## Transferability

Measurements from incompatible contexts are never pooled. Shared molecules are compared across context pairs only when at least ten are available. Rank correlation, median absolute pActivity shift, and threshold contradiction rate are reported alongside assay-type, BAO-format, and organism agreement. Organism transfer cannot be estimated in this release because the selected target contexts are all human and the activity records do not provide a consistently populated experimental-system organism field.

## Reproduction

The exact executed command is saved as `experiments/multi_target_benchmark/command.txt`. `manifest.json` records SHA-256 hashes for every normalized activity file, natural-product crossmatch, target specification, combined SQLite index, and result JSON. Repeated runs with the same inputs and seed must reproduce cohort assignments, CSVs, metrics, and plots; only generation timestamps and path-bearing manifests may differ.

## Interpretation limits

Five targets remain a small target-level sample. Publication IDs are imperfect proxies for independent laboratories. Assay descriptions do not fully encode protocol, batch, construct, or replicate differences. Natural-product evaluation cohorts remain small for some targets. Evidence scores and reliability values are heuristics, not calibrated probabilities.
