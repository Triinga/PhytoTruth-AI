import { useMemo, useState } from 'react'
import { AlertTriangle, Beaker, BrainCircuit, CheckCircle2, FlaskConical, ShieldCheck, Target } from 'lucide-react'
import { CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis } from 'recharts'
import rawResults from './data/evidence-intelligence-results.json'
import type { EvidenceIntelligenceResults } from './types'

const results = rawResults as unknown as EvidenceIntelligenceResults

function pct(value: number) {
  return `${(value * 100).toFixed(0)}%`
}

export function StudyPlannerWorkspace() {
  const [selectedId, setSelectedId] = useState(results.recommendations[0]?.molecule_id ?? '')
  const selected = results.recommendations.find(item => item.molecule_id === selectedId) ?? results.recommendations[0]
  const curve = useMemo(() => results.selective_prediction.map(row => ({
    retained: `${Math.round(row.retained_fraction * 100)}%`,
    mae: row.mae,
    coverage: row.coverage_90,
  })), [])
  const d = results.disagreement_evaluation
  const distribution = results.forward_temporal_distribution_metrics

  return <>
    <header className="page-title planner-title"><div><span className="section-kicker">STUDY-AWARE EVIDENCE INTELLIGENCE</span><h1>Don’t only predict the molecule. Plan the evidence.</h1><p>PhytoTruth estimates a future-study distribution, quantifies ten-fold disagreement risk, abstains when uncertainty is high, and proposes the independent experiment expected to reduce uncertainty most.</p></div><div className="planner-status"><BrainCircuit/><div><strong>Research-planning guardrail</strong><span>Retrospective · not autonomous</span></div></div></header>

    <section className="planner-metrics">
      <div><Target/><span>ROBUST METHOD</span><strong>Student‑t</strong><small>#1 study-grouped · #2 temporal NLL</small></div>
      <div><ShieldCheck/><span>90% COVERAGE</span><strong>{pct(distribution.coverage_90)}</strong><small>{distribution.studies} forward-time studies</small></div>
      <div><FlaskConical/><span>DISAGREEMENT AUROC</span><strong>{d.auroc.toFixed(3)}</strong><small>95% study-cluster CI {d.auroc_study_cluster_ci95[0].toFixed(3)}–{d.auroc_study_cluster_ci95[1].toFixed(3)}</small></div>
      <div><AlertTriangle/><span>POSITIVE EVENTS</span><strong>{d.disagreements}/{d.n}</strong><small>limited power · no superiority claim</small></div>
    </section>

    <div className="planner-grid">
      <section className="panel planner-panel">
        <div className="panel-head"><div><span className="section-kicker">SELECTIVE PREDICTION</span><h2>Accuracy after uncertainty-based abstention</h2></div><span className="artifact-badge">REAL FORWARD-TIME RESULTS</span></div>
        <div className="planner-chart"><ResponsiveContainer width="100%" height="100%"><LineChart data={curve} margin={{ top: 18, right: 28, left: 0, bottom: 8 }}><CartesianGrid strokeDasharray="3 5" stroke="#1b3430"/><XAxis dataKey="retained" tick={{ fill: '#78928b', fontSize: 12 }} label={{ value: 'Predictions retained', position: 'insideBottom', offset: -4, fill: '#78928b', fontSize: 12 }}/><YAxis yAxisId="mae" domain={[0, .6]} tick={{ fill: '#78928b', fontSize: 12 }}/><YAxis yAxisId="coverage" orientation="right" domain={[0, 1]} tickFormatter={pct} tick={{ fill: '#78928b', fontSize: 12 }}/><Tooltip contentStyle={{ background: '#0a1815', border: '1px solid #28473f', borderRadius: 8 }}/><Line yAxisId="mae" type="monotone" dataKey="mae" name="MAE" stroke="#6bf5bb" strokeWidth={3}/><Line yAxisId="coverage" type="monotone" dataKey="coverage" name="90% coverage" stroke="#f0ba69" strokeWidth={2}/></LineChart></ResponsiveContainer></div>
        <div className="planner-finding"><CheckCircle2/><p>Keeping the 60% least-uncertain cases reduced MAE from <strong>{curve[0].mae.toFixed(3)}</strong> to <strong>{curve[2].mae.toFixed(3)}</strong>, while observed 90% interval coverage rose from <strong>{pct(curve[0].coverage)}</strong> to <strong>{pct(curve[2].coverage)}</strong>. This is retrospective and needs independent replication.</p></div>
      </section>

      <aside className="panel planner-panel planner-inference">
        <span className="section-kicker">DISAGREEMENT CALIBRATION</span><h2>Will the next study shift ≥10×?</h2>
        <div className="risk-dial"><strong>{d.auroc.toFixed(3)}</strong><span>AUROC</span></div>
        <dl><div><dt>Brier score</dt><dd>{d.brier.toFixed(3)}</dd></div><div><dt>Average precision</dt><dd>{d.average_precision.toFixed(3)}</dd></div><div><dt>Observed prevalence</dt><dd>{pct(d.prevalence)}</dd></div><div><dt>Independent study clusters</dt><dd>{d.studies}</dd></div></dl>
        <div className="integrity-mini"><AlertTriangle/><p>The AUROC interval crosses 0.5. Three disagreement events are insufficient for a validated classifier; the probability remains a research signal.</p></div>
      </aside>
    </div>

    <section className="panel planner-recommendations">
      <div className="panel-head"><div><span className="section-kicker">EXPECTED INFORMATION GAIN</span><h2>Which experiment should be run next?</h2></div><span className="artifact-badge">OUTCOME-BLIND PRIORITY</span></div>
      <div className="recommendation-layout">
        <div className="recommendation-list">{results.recommendations.map((item, index) => <button key={item.molecule_id} className={item.molecule_id === selectedId ? 'active' : ''} onClick={() => setSelectedId(item.molecule_id)}><b>{String(index + 1).padStart(2, '0')}</b><div><strong>{item.molecule_id}</strong><span>{item.recommended_experiment}</span></div><em>{item.expected_information_gain_nats.toFixed(3)} nat</em></button>)}</div>
        {selected && <article className="experiment-prescription"><div className="prescription-head"><Beaker/><div><span>RECOMMENDED NEXT STUDY</span><h3>{selected.recommended_experiment}</h3></div></div><div className="prediction-band"><span style={{ left: `${Math.max(2, Math.min(98, (selected.predicted_future_pactivity - 4) / 6 * 100))}%` }}/><i>{selected.prediction_interval_90[0].toFixed(2)}</i><strong>predicted {selected.predicted_future_pactivity.toFixed(2)}</strong><i>{selected.prediction_interval_90[1].toFixed(2)}</i></div><div className="prescription-grid"><div><span>Disagreement probability</span><strong>{pct(selected.disagreement_probability)}</strong></div><div><span>Prior independent studies</span><strong>{selected.prior_studies}</strong></div><div><span>Protocol missingness</span><strong>{pct(selected.protocol_missing_fraction)}</strong></div><div><span>Information gain</span><strong>{selected.expected_information_gain_nats.toFixed(3)} nat</strong></div></div><p>{selected.rationale}</p><small>Pre-register the protocol and use a new publication/laboratory batch. This recommendation did not use the held-out outcome.</small></article>}
      </div>
    </section>

    <section className="experiment-conclusion"><ShieldCheck/><div><span>SCIENTIFIC BOUNDARY</span><strong>{results.scientific_conclusion}</strong><p>{results.integrity.warning}</p></div></section>
  </>
}
