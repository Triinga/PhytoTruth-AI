import { useMemo, useState } from 'react'
import { AlertTriangle, BarChart3, BookOpen, FlaskConical, Info } from 'lucide-react'
import { Bar, BarChart, CartesianGrid, Cell, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from 'recharts'
import rawResults from './data/robust-multitarget-results.json'
import type { RobustEvaluationResults, RobustMetricKey } from './types'

const results = rawResults as unknown as RobustEvaluationResults

const metrics: { key: RobustMetricKey; label: string; short: string }[] = [
  { key: 'ndcg_at_k', label: 'nDCG@10', short: 'Top-10 ranking quality' },
  { key: 'precision_at_k', label: 'Precision@10', short: 'Top-10 active fraction' },
  { key: 'enrichment_factor_at_k', label: 'EF@10', short: 'Active enrichment' },
  { key: 'kendall_tau', label: 'Kendall τ', short: 'Whole-list agreement' },
]

const ablationLabels: Record<string, string> = {
  without_support: 'Remove support',
  without_consistency: 'Remove consistency',
  without_reliability: 'Remove reliability',
  without_recent_drift: 'Remove drift',
}

function format(value: number, digits = 3) {
  return value.toFixed(digits)
}

function signed(value: number, digits = 3) {
  return `${value >= 0 ? '+' : ''}${format(value, digits)}`
}

export function RobustExperimentWorkspace() {
  const [metricKey, setMetricKey] = useState<RobustMetricKey>('ndcg_at_k')
  const comparison = results.evaluation_inference[metricKey]
  const currentMetric = metrics.find(metric => metric.key === metricKey)!
  const supportsImprovement = comparison.difference_ci95[0] > 0 && comparison.holm_adjusted_p < 0.05

  const comparisonData = useMemo(() => metrics.map(metric => ({
    metric: metric.label,
    'Median potency': results.evaluation_inference[metric.key].baseline,
    'Optimized ranker': results.evaluation_inference[metric.key].candidate,
  })), [])

  const targetData = useMemo(() => Object.entries(comparison.target_differences).map(([target, difference]) => ({
    target: target.replace('CHEMBL', ''), difference,
  })), [comparison])

  const ablationData = useMemo(() => Object.entries(results.ablations).flatMap(([name, ablation]) => {
    if (!ablation.full_minus_ablation) return []
    return [{ variant: ablationLabels[name] ?? name, value: ablation.full_minus_ablation[metricKey].difference }]
  }), [metricKey])

  const currentHeuristic = results.current_heuristic_evaluation_inference.ndcg_at_k

  return <>
    <header className="page-title experiment-title"><div><span className="section-kicker">LOCKED MULTI-TARGET EVALUATION</span><h1>Evidence first, even when the result is negative</h1><p>Twenty eligible biological targets were split before model fitting. The latest publication year for each molecule was held out, while model and parameter selection used development targets only.</p></div><div className="experiment-verdict"><AlertTriangle/><div><strong>Superiority not demonstrated</strong><span>Frozen protocol · seed {results.seed}</span></div></div></header>

    <section className="experiment-context-strip">
      <div><span>ELIGIBLE TARGETS</span><strong>{results.coverage.eligible_targets}</strong><small>{results.coverage.discovered_targets} discovered locally</small></div>
      <div><span>UNTOUCHED EVALUATION TARGETS</span><strong>{results.coverage.evaluation_targets}</strong><small>{results.coverage.development_targets} targets used for development</small></div>
      <div><span>HELD-OUT MOLECULES</span><strong>{results.coverage.evaluation_molecules.toLocaleString()}</strong><small>latest-publication evidence only</small></div>
      <div><span>TARGET BOOTSTRAPS</span><strong>{results.bootstrap_iterations.toLocaleString()}</strong><small>independence unit: biological target</small></div>
    </section>

    <div className="experiment-controls" aria-label="Experiment filters">
      <div><label htmlFor="metric-select">Metric focus</label><select id="metric-select" value={metricKey} onChange={event => setMetricKey(event.target.value as RobustMetricKey)}>{metrics.map(item => <option value={item.key} key={item.key}>{item.label} · {item.short}</option>)}</select><small>Top-K metrics use k = {results.k}</small></div>
      <div><label>DEVELOPMENT-ONLY SELECTION</label><strong className="control-value">{results.development_selection.candidate_config}</strong><small>{results.development_selection.passed_development_gate ? 'Passed' : 'Did not pass'} the frozen two-metric gate</small></div>
      <div className="cohort-size"><span>DEPLOYED RANKING</span><strong>Median</strong><small>guardrail retained after negative evaluation</small></div>
    </div>

    <div className="experiment-grid">
      <section className="panel experiment-chart-panel">
        <div className="panel-head"><div><span className="section-kicker">TARGET-MACRO PERFORMANCE</span><h2>Held-out ranking performance</h2></div><BarChart3 size={19}/></div>
        <div className="experiment-chart"><ResponsiveContainer width="100%" height="100%"><BarChart data={comparisonData} margin={{ top: 14, right: 12, bottom: 8, left: -14 }}><CartesianGrid strokeDasharray="3 5" stroke="#1b3430" vertical={false}/><XAxis dataKey="metric" tick={{ fill: '#78928b', fontSize: 12 }} axisLine={{ stroke: '#29443e' }}/><YAxis tick={{ fill: '#78928b', fontSize: 12 }} axisLine={false}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8, fontSize: 12 }}/><Legend wrapperStyle={{ fontSize: 12 }}/><Bar dataKey="Median potency" fill="#f0ba69" radius={[4, 4, 0, 0]}/><Bar dataKey="Optimized ranker" fill="#6bf5bb" radius={[4, 4, 0, 0]}/></BarChart></ResponsiveContainer></div>
      </section>

      <aside className="panel metric-inspector">
        <span className="section-kicker">{currentMetric.label} · TARGET-LEVEL INFERENCE</span>
        <div className="metric-versus"><div><span>Median potency</span><strong>{format(comparison.baseline)}</strong></div><i>vs</i><div><span>Optimized</span><strong>{format(comparison.candidate)}</strong></div></div>
        <div className={`difference-callout ${supportsImprovement ? 'significant' : 'inconclusive'}`}><span>MACRO DIFFERENCE</span><strong>{signed(comparison.difference)}</strong><small>95% target-bootstrap CI {signed(comparison.difference_ci95[0])} to {signed(comparison.difference_ci95[1])}</small></div>
        <div className="significance-row"><Info/><div><strong>{supportsImprovement ? 'Improvement supported' : 'No reliable improvement'}</strong><span>exact p = {comparison.target_sign_flip_p.toFixed(6)} · Holm p = {comparison.holm_adjusted_p.toFixed(6)}</span></div></div>
        <p>Paired Hedges g = {format(comparison.hedges_g_paired)}. Target wins/ties/losses: {comparison.target_wins}/{comparison.target_ties}/{comparison.target_losses}. Molecules are not treated as independent biological replications.</p>
      </aside>
    </div>

    <section className="panel ablation-panel">
      <div className="panel-head"><div><span className="section-kicker">TARGET-BY-TARGET EFFECT</span><h2>{currentMetric.label}: optimized minus median</h2></div><span className="artifact-badge">ACTUAL ARTIFACT</span></div>
      <div className="ablation-layout"><div className="ablation-chart"><ResponsiveContainer width="100%" height="100%"><BarChart data={targetData} margin={{ top: 8, right: 20, bottom: 8, left: 2 }}><CartesianGrid strokeDasharray="3 5" stroke="#1b3430" vertical={false}/><XAxis dataKey="target" tick={{ fill: '#8da29d', fontSize: 12 }} axisLine={false}/><YAxis tick={{ fill: '#78928b', fontSize: 12 }}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8, fontSize: 12 }}/><Bar dataKey="difference" radius={[4, 4, 0, 0]}>{targetData.map(item => <Cell key={item.target} fill={item.difference >= 0 ? '#6bf5bb' : '#ff7185'}/>)}</Bar></BarChart></ResponsiveContainer></div><div className="negative-result"><FlaskConical/><div><strong>Power boundary, reported before claims</strong><p>With {results.power_analysis.independent_evaluation_targets} evaluation targets, the smallest attainable two-sided exact p-value is {results.power_analysis.minimum_attainable_two_sided_exact_p.toFixed(6)}. Therefore p ≤ 0.001 is impossible in this locked holdout, regardless of effect size.</p></div></div></div>
    </section>

    <section className="panel ablation-panel">
      <div className="panel-head"><div><span className="section-kicker">FROZEN ABLATIONS</span><h2>Full optimized model minus component-removed model</h2></div><span className="artifact-badge">SECONDARY ANALYSIS</span></div>
      <div className="ablation-layout"><div className="ablation-chart"><ResponsiveContainer width="100%" height="100%"><BarChart data={ablationData} layout="vertical" margin={{ top: 8, right: 24, bottom: 8, left: 30 }}><CartesianGrid strokeDasharray="3 5" stroke="#1b3430" horizontal={false}/><XAxis type="number" tick={{ fill: '#78928b', fontSize: 12 }}/><YAxis dataKey="variant" type="category" width={112} tick={{ fill: '#8da29d', fontSize: 12 }} axisLine={false}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8, fontSize: 12 }}/><Bar dataKey="value" fill="#6bf5bb" radius={[0, 4, 4, 0]}/></BarChart></ResponsiveContainer></div><div className="negative-result"><FlaskConical/><div><strong>No component earns a causal claim</strong><p>Positive bars favor the full model; negative bars favor removal. All exact tests are target-level and multiplicity-adjusted. Feature groups absent from the selected model are marked not applicable in the artifact.</p></div></div></div>
    </section>

    <section className="experiment-conclusion"><BookOpen/><div><span>SCIENTIFIC CONCLUSION</span><strong>{results.scientific_conclusion.summary}</strong><p>The prior heuristic changed nDCG@10 by {signed(currentHeuristic.difference)} (95% CI {signed(currentHeuristic.difference_ci95[0])} to {signed(currentHeuristic.difference_ci95[1])}). This is a retrospective reused-database benchmark, not prospective independent-laboratory confirmation.</p></div></section>
  </>
}
