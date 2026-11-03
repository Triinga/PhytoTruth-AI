import { useMemo, useState } from 'react'
import { AlertTriangle, BarChart3, BookOpen, CheckCircle2, FlaskConical, Info } from 'lucide-react'
import { Bar, BarChart, CartesianGrid, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from 'recharts'
import rawResults from './data/evidence-reversal-results.json'
import type { ReversalResults } from './types'

export { RobustExperimentWorkspace as ExperimentWorkspace } from './RobustExperimentWorkspace'

const results = rawResults as unknown as ReversalResults
type CohortKey = keyof ReversalResults['cohorts']
type MetricKey = 'ndcg_at_k' | 'precision_at_k' | 'enrichment_factor_at_k' | 'kendall_tau'

const cohorts: { key: CohortKey; label: string; note: string }[] = [
  { key: 'latest_all', label: 'Latest-year · all', note: 'Primary molecule-specific temporal holdout' },
  { key: 'latest_natural', label: 'Latest-year · NPASS', note: 'Verified natural-product subset' },
  { key: 'global_all', label: 'Global 2015 · all', note: 'Strict shared publication-time cutoff' },
  { key: 'global_natural', label: 'Global 2015 · NPASS', note: 'Strict natural-product sensitivity analysis' },
]

const metrics: { key: MetricKey; label: string; short: string }[] = [
  { key: 'ndcg_at_k', label: 'nDCG@10', short: 'Ranking quality' },
  { key: 'precision_at_k', label: 'Precision@10', short: 'Top-10 active fraction' },
  { key: 'enrichment_factor_at_k', label: 'EF@10', short: 'Active enrichment' },
  { key: 'kendall_tau', label: 'Kendall τ', short: 'Global ordering' },
]

function format(value: number) {
  return value.toFixed(3)
}

export function LegacyExperimentWorkspace() {
  const [cohortKey, setCohortKey] = useState<CohortKey>('latest_all')
  const [metricKey, setMetricKey] = useState<MetricKey>('ndcg_at_k')
  const cohort = results.cohorts[cohortKey]
  const comparison = cohort.baseline_vs_phytotruth[metricKey]
  const currentMetric = metrics.find(metric => metric.key === metricKey)!
  const significant = comparison.difference_ci95[0] > 0 || comparison.difference_ci95[1] < 0

  const comparisonData = useMemo(() => metrics.map(metric => ({
    metric: metric.label,
    Baseline: cohort.baseline_vs_phytotruth[metric.key].baseline,
    PhytoTruth: cohort.baseline_vs_phytotruth[metric.key].candidate,
  })), [cohort])

  const ablationData = useMemo(() => {
    const selected = cohort.evaluation.metrics
    return [
      ['Potency max', selected.baseline_max_potency?.[metricKey]],
      ['Potency median', selected.median_potency?.[metricKey]],
      ['Full heuristic', selected.phytotruth_full?.[metricKey]],
      ['No reliability', selected.ablate_reliability?.[metricKey]],
      ['No support', selected.ablate_support?.[metricKey]],
      ['No consistency', selected.ablate_explicit_consistency?.[metricKey]],
    ].filter((row): row is [string, number] => typeof row[1] === 'number').map(([variant, value]) => ({ variant, value }))
  }, [cohort, metricKey])

  return <>
    <header className="page-title experiment-title"><div><span className="section-kicker">EXECUTED EVIDENCE REVERSAL EXPERIMENT</span><h1>The result is informative—even when it is negative</h1><p>Explore the actual leakage-resistant CHEMBL335 benchmark. Rankings were frozen from earlier evidence and scored only against molecule-held-out future evidence.</p></div><div className="experiment-verdict"><AlertTriangle/><div><strong>Advantage not demonstrated</strong><span>Prespecified analysis · seed {results.seed}</span></div></div></header>

    <section className="experiment-context-strip">
      <div><span>COMPATIBLE MEASUREMENTS</span><strong>{results.selected_measurements.toLocaleString()}</strong></div>
      <div><span>INDEPENDENT MOLECULES</span><strong>{results.selected_molecules.toLocaleString()}</strong></div>
      <div><span>NPASS-VERIFIED</span><strong>{results.verified_natural_product_molecules_in_context.toLocaleString()}</strong></div>
      <div><span>RESAMPLES</span><strong>{results.bootstrap_iterations.toLocaleString()}</strong><small>bootstrap + {results.permutation_iterations.toLocaleString()} permutations</small></div>
    </section>

    <div className="experiment-controls" aria-label="Experiment filters">
      <div><label htmlFor="cohort-select">Evaluation cohort</label><select id="cohort-select" value={cohortKey} onChange={event => setCohortKey(event.target.value as CohortKey)}>{cohorts.map(item => <option value={item.key} key={item.key}>{item.label}</option>)}</select><small>{cohorts.find(item => item.key === cohortKey)?.note}</small></div>
      <div><label htmlFor="metric-select">Metric focus</label><select id="metric-select" value={metricKey} onChange={event => setMetricKey(event.target.value as MetricKey)}>{metrics.map(item => <option value={item.key} key={item.key}>{item.label} · {item.short}</option>)}</select><small>Top-K uses k = {cohort.evaluation.k}</small></div>
      <div className="cohort-size"><span>HELD-OUT COHORT</span><strong>{cohort.evaluation.n_molecules}</strong><small>{cohort.evaluation.heldout_active_molecules} active at pChEMBL ≥ {cohort.evaluation.active_threshold}</small></div>
    </div>

    <div className="experiment-grid">
      <section className="panel experiment-chart-panel">
        <div className="panel-head"><div><span className="section-kicker">BASELINE VS PHYTOTRUTH</span><h2>Held-out ranking performance</h2></div><BarChart3 size={19}/></div>
        <div className="experiment-chart"><ResponsiveContainer width="100%" height="100%"><BarChart data={comparisonData} margin={{ top: 14, right: 12, bottom: 8, left: -14 }}>
          <CartesianGrid strokeDasharray="3 5" stroke="#1b3430" vertical={false}/><XAxis dataKey="metric" tick={{ fill: '#78928b', fontSize: 12 }} axisLine={{ stroke: '#29443e' }}/><YAxis tick={{ fill: '#78928b', fontSize: 12 }} axisLine={false}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8, fontSize: 12 }}/><Legend wrapperStyle={{ fontSize: 12 }}/><Bar dataKey="Baseline" fill="#f0ba69" radius={[4, 4, 0, 0]}/><Bar dataKey="PhytoTruth" fill="#6bf5bb" radius={[4, 4, 0, 0]}/>
        </BarChart></ResponsiveContainer></div>
      </section>

      <aside className="panel metric-inspector">
        <span className="section-kicker">{currentMetric.label} INFERENCE</span>
        <div className="metric-versus"><div><span>Potency-only</span><strong>{format(comparison.baseline)}</strong></div><i>vs</i><div><span>PhytoTruth</span><strong>{format(comparison.candidate)}</strong></div></div>
        <div className={`difference-callout ${significant ? 'significant' : 'inconclusive'}`}><span>PAIRED DIFFERENCE</span><strong>{comparison.difference >= 0 ? '+' : ''}{format(comparison.difference)}</strong><small>95% CI {format(comparison.difference_ci95[0])} to {format(comparison.difference_ci95[1])}</small></div>
        <div className="significance-row">{significant ? <CheckCircle2/> : <Info/>}<div><strong>{significant ? 'Interval excludes zero' : 'No significant advantage'}</strong><span>paired randomization p = {comparison.paired_randomization_p.toFixed(4)}</span></div></div>
        <p>Confidence intervals are molecule-bootstrap intervals. A higher point estimate does not demonstrate improvement when its paired interval includes zero.</p>
      </aside>
    </div>

    <section className="panel ablation-panel">
      <div className="panel-head"><div><span className="section-kicker">ABLATION STUDY</span><h2>Do the reliability components add value?</h2></div><span className="artifact-badge">ACTUAL ARTIFACT</span></div>
      <div className="ablation-layout"><div className="ablation-chart"><ResponsiveContainer width="100%" height="100%"><BarChart data={ablationData} layout="vertical" margin={{ top: 8, right: 24, bottom: 8, left: 28 }}><CartesianGrid strokeDasharray="3 5" stroke="#1b3430" horizontal={false}/><XAxis type="number" tick={{ fill: '#78928b', fontSize: 12 }}/><YAxis dataKey="variant" type="category" width={96} tick={{ fill: '#8da29d', fontSize: 12 }} axisLine={false}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8, fontSize: 12 }}/><Bar dataKey="value" fill="#6bf5bb" radius={[0, 4, 4, 0]}/></BarChart></ResponsiveContainer></div><div className="negative-result"><FlaskConical/><div><strong>Negative finding, reported plainly</strong><p>{results.scientific_conclusion.ablation_finding}</p></div></div></div>
    </section>

    <section className="experiment-conclusion"><BookOpen/><div><span>SCIENTIFIC CONCLUSION</span><strong>{results.scientific_conclusion.summary}</strong><p>{results.scientific_conclusion.natural_product_finding}</p></div></section>
  </>
}
