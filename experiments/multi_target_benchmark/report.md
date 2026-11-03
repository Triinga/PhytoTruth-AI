# PhytoTruth multi-target evidence benchmark

Generated 2026-10-20T15:17:54.992129+00:00 from real normalized ChEMBL IC50 measurements and structure-verified NPASS/COCONUT membership.

## Design

Five human single-protein targets were prespecified before outcomes were inspected. A target was included only if its largest compatible target/organism/endpoint/assay-type/BAO context had at least 30 molecules with strictly earlier, document-disjoint training publications and at least 5 structure-verified natural products. Target selection used coverage, never activity values or ranking results.

Each included target used a latest-publication temporal holdout followed by a deterministic molecule-ID development/evaluation split stratified by natural-product membership. Methods were selected only from development results. Evaluation results were opened after the development gate was frozen. The prior median-potency guardrail remained the default unless an evidence-aware method improved both macro nDCG and Kendall tau-b with intervals excluding zero and p<0.05.

## Coverage

- CHEMBL220 (Acetylcholinesterase): 105 evaluation molecules (17 natural products); median/current/lower-CI nDCG@10 = 0.947/0.944/0.948.
- CHEMBL203 (Epidermal growth factor receptor): 326 evaluation molecules (12 natural products); median/current/lower-CI nDCG@10 = 0.959/0.927/0.906.
- CHEMBL230 (Prostaglandin G/H synthase 2): 194 evaluation molecules (12 natural products); median/current/lower-CI nDCG@10 = 0.976/0.976/0.976.
- CHEMBL284 (Dipeptidyl peptidase 4): 73 evaluation molecules (6 natural products); median/current/lower-CI nDCG@10 = 0.892/0.840/0.892.
- CHEMBL335 (Tyrosine-protein phosphatase non-receptor type 1): 64 evaluation molecules (14 natural products); median/current/lower-CI nDCG@10 = 0.929/0.929/0.929.

The combined natural-product index contains 702036 unique standardized structures with 201971 NPASS and 695128 COCONUT provenance records. Cross-source duplicate structures retain both provenance records.

## Independent evaluation

The development gate selected `median_potency`. Across targets, current PhytoTruth minus median nDCG was -0.017 (95% CI -0.043 to 0.001; p=0.0569), and Kendall difference was -0.024. Lower-confidence-bound minus median nDCG was -0.010 (95% CI -0.048 to 0.002; p=0.0480). These results do not support superiority unless the recorded development gate and independent evaluation both do so.

## Transferability

`assay_transferability.csv` compares shared compounds across incompatible contexts without pooling their measurements. It reports rank correlation, median absolute pActivity shift, and activity-threshold contradictions, stratified by organism, assay type, and BAO format agreement. These comparisons are descriptive and do not assume that evidence transfers across contexts.

Organism transfer was not estimable: the prespecified targets and all retained target contexts were human. ChEMBL activity records used here expose target organism but not a consistently populated experimental-system organism field, so no cross-organism claim is made.

## Limitations

- Evaluation is literature-temporal and molecule-partitioned, but not a prospective or external laboratory validation.
- Publication IDs do not prove laboratory independence; protocols, batches, and replicate structure remain incompletely recorded.
- Natural-product subgroup power varies by target despite combining NPASS and COCONUT.
- Target-level inference has only 5 biological targets, so hierarchical intervals remain imprecise.
- Held-out bioactivity is not evidence of clinical efficacy or safety.
