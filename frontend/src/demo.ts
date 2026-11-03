import type { Dashboard } from './types'

// This fallback is a checked-in rendering of backend/data/ptp1b_snapshot.json.
// Values are derived from ChEMBL API records retrieved 2026-10-20, never invented.
export const fallbackDashboard: Dashboard = {
  metadata: {
    source: 'ChEMBL REST API', target_chembl_id: 'CHEMBL335',
    target_name: 'Tyrosine-protein phosphatase non-receptor type 1 (PTP1B)',
    organism: 'Homo sapiens', standard_type: 'IC50', retrieved_at: '2026-10-20',
    query: 'target_chembl_id=CHEMBL335&standard_type=IC50&standard_units=nM&pchembl_value__isnull=false',
    note: 'Small, traceable demonstration slice. Natural-product status is shown separately from activity evidence.',
  },
  summary: { molecules: 5, measurements: 17, independent_assays: 11, flagged_candidates: 2 },
  candidates: [
    { id:'CHEMBL430266', name:'Calyculin A', smiles:'COC[C@@H]([C@H](O)[C@H](O)C(=O)NCC[C@H](C)c1nc(/C=C/C[C@@H]2O[C@]3(C[C@@H](O)[C@@H]2C)O[C@H]([C@H](C[C@H](O)[C@H](C)[C@H](O)[C@H](C)/C=C(C)/C(C)=C/C=C/C(C)=C\\C#N)OC)[C@H](OP(=O)(O)O)C3(C)C)co1)N(C)C',natural_product_status:'known natural product; verify organism provenance through COCONUT/NPASS import',measurements:3,independent_assays:3,documents:3,median_pchembl:9.3,max_pchembl:9.52,mad:.22,pchembl_range:.6,contradiction_rate:0,duplicate_rate:0,organism_count:1,reliability:.626,evidence_score:9.335,confidence_interval:[8.92,9.52],flags:[],potency_rank:1,evidence_rank:1,rank_change:0 },
    { id:'CHEMBL280487', name:'Okadaic acid', smiles:'C=C1[C@@H](O)[C@@H]2O[C@]3(CC[C@H](/C=C/[C@@H](C)[C@@H]4CC(C)=C[C@@]5(O[C@H](C[C@@](C)(O)C(=O)O)CC[C@H]5O)O4)O3)CC[C@H]2O[C@@H]1[C@@H](O)C[C@H](C)[C@H]1O[C@@]2(CCCCO2)CC[C@H]1C',natural_product_status:'known natural product; verify organism provenance through COCONUT/NPASS import',measurements:4,independent_assays:4,documents:3,median_pchembl:7.845,max_pchembl:8.91,mad:.845,pchembl_range:3.02,contradiction_rate:.5,duplicate_rate:0,organism_count:1,reliability:.288,evidence_score:6.823,confidence_interval:[5.89,8.91],flags:['high assay spread','activity threshold reversal'],potency_rank:2,evidence_rank:2,rank_change:0 },
    { id:'CHEMBL96952', name:'CHEMBL96952', smiles:'Nc1nc(N)c2nc(CN3CCN(Cc4ccccc4)CC3)nnc2n1',natural_product_status:'not yet cross-matched',measurements:3,independent_assays:2,documents:1,median_pchembl:5.1,max_pchembl:5.21,mad:.11,pchembl_range:1.2,contradiction_rate:0,duplicate_rate:0,organism_count:1,reliability:.499,evidence_score:5.082,confidence_interval:[4.01,5.21],flags:[],potency_rank:4,evidence_rank:3,rank_change:1 },
    { id:'CHEMBL420456', name:'CHEMBL420456', smiles:'Nc1nc(N)c2nc(CN3CCN(Cc4ccc(-c5ccccc5)cc4)CC3)nnc2n1',natural_product_status:'not yet cross-matched',measurements:3,independent_assays:2,documents:1,median_pchembl:5.1,max_pchembl:5.54,mad:.35,pchembl_range:.79,contradiction_rate:0,duplicate_rate:0,organism_count:1,reliability:.402,evidence_score:4.975,confidence_interval:[4.75,5.54],flags:[],potency_rank:3,evidence_rank:4,rank_change:-1 },
    { id:'CHEMBL89067', name:'CHEMBL89067', smiles:'CSc1ccnc(C[S+]([O-])c2nc3cc(OC(C)C)ccc3[nH]2)c1',natural_product_status:'not yet cross-matched',measurements:4,independent_assays:4,documents:2,median_pchembl:4.26,max_pchembl:4.74,mad:0,pchembl_range:.48,contradiction_rate:0,duplicate_rate:.25,organism_count:1,reliability:.731,evidence_score:4.492,confidence_interval:[4.26,4.74],flags:['potential duplicate'],potency_rank:5,evidence_rank:5,rank_change:0 },
  ],
}

