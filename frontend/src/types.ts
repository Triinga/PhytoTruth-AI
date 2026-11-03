export type Candidate = {
  id: string
  name: string
  smiles: string | null
  natural_product_status: string
  measurements: number
  independent_assays: number
  documents: number
  median_pchembl: number
  max_pchembl: number
  mad: number
  pchembl_range: number
  contradiction_rate: number
  duplicate_rate: number
  organism_count: number
  reliability: number
  evidence_score: number
  confidence_interval: [number, number]
  flags: string[]
  potency_rank: number
  evidence_rank: number
  rank_change: number
}

export type Dashboard = {
  metadata: {
    source: string
    target_chembl_id: string
    target_name: string
    organism: string
    standard_type: string
    retrieved_at: string
    query: string
    note: string
  }
  summary: { molecules: number; measurements: number; independent_assays: number; flagged_candidates: number }
  candidates: Candidate[]
}

export type PredictionInterval = {
  lower: number
  upper: number
  confidence_level: number
  method: string
}

export type ApplicabilityDomain = {
  inside: boolean
  status: 'IN_DOMAIN_RESEARCH_ONLY' | 'OUTSIDE_DOMAIN_DO_NOT_RELY'
  nearest_training_tanimoto: number
  threshold: number
  warning: string | null
}

export type MolecularDescriptors = {
  molecular_weight: number
  logp: number
  tpsa: number
  hbd: number
  hba: number
  rotatable_bonds: number
  fraction_csp3: number
  ring_count: number
  heavy_atoms: number
  formal_charge: number
}

export type PredictionResponse = {
  canonical_smiles: string
  predicted_pchembl: number
  prediction_interval: PredictionInterval
  applicability_domain: ApplicabilityDomain
  descriptors: MolecularDescriptors
  model: {
    name: string
    version: string
    format: string
    target: {
      assay_type: string | null
      bao_format: string | null
      standard_type: string | null
      target_chembl_id: string | null
      target_organism: string | null
    }
    training_structures: number
    training_scaffolds: number
  }
  kind: 'predicted'
  limitations: string[]
}

export type ScreeningTargetPrediction = PredictionResponse & {
  target_chembl_id: string
  target_name: string
  weight: number
  direction: string
  target_score: number
}

export type ScreeningArtifactCandidate = {
  rank: number
  standard_inchi_key: string
  canonical_smiles: string
  names: string | null
  reported_source_organisms: string | null
  natural_product_provenance: string | null
  disease_id: string
  disease_name: string
  prioritization_score: number
  all_targets_inside_domain: boolean
  structural_alerts: string[]
  descriptors: MolecularDescriptors
  target_predictions: ScreeningTargetPrediction[]
  research_status: 'PRIORITIZED_FOR_EXPERIMENTAL_VALIDATION' | 'OUT_OF_DOMAIN_REVIEW_REQUIRED'
}

export type ScreeningArtifactResponse = {
  candidates: ScreeningArtifactCandidate[]
  count: number
  metadata: null | {
    created_at: string
    screened_compounds: number
    property_filtered_compounds: number
    prediction_failed_compounds: number
    ranked_candidates: number
    top_k: number
    out_of_domain_quota: number
    in_domain_candidates: number
    out_of_domain_candidates: number
    source_filter: string | null
    max_compounds: number | null
    disease: { id?: string; name?: string }
    targets: Array<{ target_chembl_id: string; target_name?: string; direction?: string; weight?: number }>
    filters: Record<string, number>
    score_semantics: string
    limitations: string[]
    specifications_sha256: string
    natural_product_database_sha256: string
  }
  semantics: string
}

export type MetricComparison = {
  baseline: number
  baseline_ci95: [number, number]
  candidate: number
  candidate_ci95: [number, number]
  difference: number
  difference_ci95: [number, number]
  paired_randomization_p: number
}

export type ReversalCohort = {
  baseline_vs_phytotruth: Record<'ndcg_at_k' | 'precision_at_k' | 'enrichment_factor_at_k' | 'kendall_tau', MetricComparison>
  evaluation: {
    active_threshold: number
    heldout_active_molecules: number
    k: number
    n_molecules: number
    metrics: Record<string, Record<'ndcg_at_k' | 'precision_at_k' | 'enrichment_factor_at_k' | 'kendall_tau', number>>
  }
}

export type ReversalResults = {
  generated_at: string
  seed: number
  selected_measurements: number
  selected_molecules: number
  verified_natural_product_molecules_in_context: number
  bootstrap_iterations: number
  permutation_iterations: number
  cohorts: Record<'global_all' | 'global_natural' | 'latest_all' | 'latest_natural', ReversalCohort>
  scientific_conclusion: {
    status: string
    summary: string
    natural_product_finding: string
    ablation_finding: string
  }
}

export type RobustMetricComparison = {
  baseline: number
  candidate: number
  difference: number
  difference_ci95: [number, number]
  target_sign_flip_p: number
  holm_adjusted_p: number
  hedges_g_paired: number
  target_wins: number
  target_ties: number
  target_losses: number
  target_differences: Record<string, number>
}

export type RobustMetricKey = 'ndcg_at_k' | 'precision_at_k' | 'enrichment_factor_at_k' | 'kendall_tau'

export type RobustEvaluationResults = {
  generated_at: string
  seed: number
  k: number
  bootstrap_iterations: number
  independence_unit: string
  coverage: {
    discovered_targets: number
    eligible_targets: number
    development_targets: number
    evaluation_targets: number
    development_molecules: number
    evaluation_molecules: number
    development_target_ids: string[]
    evaluation_target_ids: string[]
  }
  development_selection: {
    candidate_config: string
    passed_development_gate: boolean
    deployed_method: string
  }
  evaluation_inference: Record<RobustMetricKey, RobustMetricComparison>
  current_heuristic_evaluation_inference: Record<RobustMetricKey, RobustMetricComparison>
  power_analysis: {
    independent_evaluation_targets: number
    minimum_attainable_two_sided_exact_p: number
    observed_paired_standardized_effect: number
    approximate_targets_for_80_percent_power_at_alpha_0_05: number | null
    note: string
  }
  ablations: Record<string, {
    removed_features: string[]
    status: string
    full_minus_ablation: Record<RobustMetricKey, RobustMetricComparison> | null
  }>
  scientific_conclusion: {
    status: string
    superiority_claim: boolean
    deployed_ranking: string
    summary: string
    confirmatory_status: string
  }
}

export type EvidenceIntelligenceResults = {
  experiment: string
  seed: number
  source_status: string
  selection: {
    rule: string
    selected_method: string
    selected: {
      method: string
      ranks: Record<'study_grouped' | 'forward_temporal', number>
      worst_nll_rank: number
      mean_nll_rank: number
      mean_calibration_error: number
    }
  }
  forward_temporal_distribution_metrics: {
    n: number
    compounds: number
    studies: number
    mae: number
    rmse: number
    nll: number
    coverage_80: number
    coverage_90: number
    mean_width_90: number
    mean_absolute_calibration_error: number
  }
  study_grouped_distribution_metrics: {
    n: number
    studies: number
    mae: number
    nll: number
    coverage_90: number
  }
  disagreement_evaluation: {
    n: number
    studies: number
    disagreements: number
    prevalence: number
    brier: number
    brier_study_cluster_ci95: [number, number]
    auroc: number
    auroc_study_cluster_ci95: [number, number]
    average_precision: number
    average_precision_study_cluster_ci95: [number, number]
  }
  selective_prediction: Array<{
    retained_fraction: number
    retained_n: number
    abstained_n: number
    mae: number
    coverage_90: number
    mean_interval_width_90: number
  }>
  recommendations: Array<{
    target_chembl_id: string
    molecule_id: string
    prior_median_pactivity: number
    predicted_future_pactivity: number
    prediction_interval_90: [number, number]
    disagreement_probability: number
    prior_studies: number
    protocol_missing_fraction: number
    expected_information_gain_nats: number
    recommended_experiment: string
    rationale: string
    decision_status: string
  }>
  integrity: {
    heldout_outcome_used_for_recommendation_priority: boolean
    confirmatory_status: string
    deployment_status: string
    warning: string
  }
  scientific_conclusion: string
}

