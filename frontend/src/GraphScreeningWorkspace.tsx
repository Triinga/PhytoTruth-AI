import { useEffect, useMemo, useRef, useState } from 'react'
import {
  AlertTriangle, Atom, CheckCircle2, ChevronRight, CircleAlert, Database,
  Dna, ExternalLink, Filter, FlaskConical, GitBranch, LoaderCircle,
  Network, Play, Search, ShieldCheck, Sparkles, Target, TestTube2,
} from 'lucide-react'
import { ApiError, apiDocsUrl, getScreeningCandidates, predictMolecule } from './api'
import { MoleculeViewer } from './MoleculeViewer'
import type {
  Candidate, Dashboard, MolecularDescriptors, PredictionResponse,
  ScreeningArtifactCandidate, ScreeningArtifactResponse,
} from './types'

type ScreeningRow = {
  rank: number
  key: string
  id: string
  name: string
  smiles: string
  source: string
  organisms: string
  evidence: string
  measurements: number | null
  assays: number | null
  score: number
  prediction: PredictionResponse
  alerts: string[]
  descriptors: MolecularDescriptors
}

type GraphNodeKey = 'source' | 'compound' | 'target' | 'disease' | 'assay' | 'publication'

const nodeCopy: Record<GraphNodeKey, { label: string; edge: string }> = {
  source: { label: 'Natural-product source', edge: 'VERIFIED_BY_SOURCE' },
  compound: { label: 'Compound', edge: 'PREDICTED_AGAINST' },
  target: { label: 'Protein target', edge: 'ASSOCIATED_WITH_DISEASE' },
  disease: { label: 'Disease', edge: 'CURATED_RELATIONSHIP' },
  assay: { label: 'Assay context', edge: 'MEASURES' },
  publication: { label: 'Evidence records', edge: 'REPORTS' },
}

function demoScore(prediction: PredictionResponse): number {
  const width = prediction.prediction_interval.upper - prediction.prediction_interval.lower
  return prediction.predicted_pchembl
    + (prediction.applicability_domain.inside ? .35 : -1.25)
    + .25 * prediction.applicability_domain.nearest_training_tanimoto
    - .15 * Math.max(0, width)
}

function isVerifiedNaturalProduct(candidate: Candidate): boolean {
  return candidate.natural_product_status.toLowerCase().includes('known natural product')
}

function artifactRow(candidate: ScreeningArtifactCandidate): ScreeningRow | null {
  const target = candidate.target_predictions[0]
  if (!target) return null
  const provenanceLabel = candidate.natural_product_provenance?.split(',')[0]?.replace(':', ' · ')
  return {
    rank: candidate.rank,
    key: candidate.standard_inchi_key,
    id: candidate.standard_inchi_key,
    name: candidate.names?.split(',')[0] || provenanceLabel || candidate.standard_inchi_key,
    smiles: candidate.canonical_smiles,
    source: candidate.natural_product_provenance || 'NPASS / COCONUT',
    organisms: candidate.reported_source_organisms || 'Not reported in structure release',
    evidence: 'Precomputed checksum-tracked screening artifact',
    measurements: null,
    assays: null,
    score: candidate.prioritization_score,
    prediction: target,
    alerts: candidate.structural_alerts,
    descriptors: candidate.descriptors,
  }
}

function Stat({ label, value, note }: { label: string; value: string | number; note: string }) {
  return <div className="gs-stat"><span>{label}</span><strong>{value}</strong><small>{note}</small></div>
}

function GraphLegend() {
  return <div className="gs-legend">
    <span><i className="measured"/> measured</span>
    <span><i className="predicted"/> predicted</span>
    <span><i className="curated"/> curated</span>
  </div>
}

function EvidenceGraph({ row }: { row: ScreeningRow }) {
  const [active, setActive] = useState<GraphNodeKey>('compound')
  const target = row.prediction.model.target
  const sourceNodeLabel = row.source.toLowerCase().includes('known natural product')
    ? 'NPASS / COCONUT'
    : row.source.split(',')[0]
  const details: Record<GraphNodeKey, { title: string; subtitle: string; evidence: string }> = {
    source: {
      title: row.source,
      subtitle: row.organisms,
      evidence: row.evidence,
    },
    compound: {
      title: row.name,
      subtitle: row.id,
      evidence: `Canonical structure: ${row.prediction.canonical_smiles}`,
    },
    target: {
      title: 'PTP1B',
      subtitle: target.target_chembl_id || 'CHEMBL335',
      evidence: `${target.target_organism || 'Homo sapiens'} · ${target.standard_type || 'IC50'} model`,
    },
    disease: {
      title: 'Type 2 diabetes mellitus',
      subtitle: 'MONDO:0005148',
      evidence: 'Curated research-target relationship from the graph specification; not a causal or clinical claim.',
    },
    assay: {
      title: 'Biochemical IC₅₀ assay',
      subtitle: target.bao_format || 'BAO_0000357',
      evidence: `${row.assays ?? 'Multiple'} traceable assay context${row.assays === 1 ? '' : 's'}`,
    },
    publication: {
      title: `${row.measurements ?? 'Multiple'} evidence record${row.measurements === 1 ? '' : 's'}`,
      subtitle: 'ChEMBL provenance',
      evidence: 'Activity, assay, document and source identifiers remain separate from ML predictions.',
    },
  }
  const node = details[active]

  return <section className="panel gs-graph-panel">
    <div className="panel-head gs-panel-head">
      <div><span className="section-kicker">SEMANTIC EVIDENCE GRAPH</span><h2>Trace the candidate, not just the score</h2></div>
      <GraphLegend/>
    </div>
    <div className="gs-graph" aria-label="Interactive semantic knowledge graph">
      <svg viewBox="0 0 1000 410" preserveAspectRatio="none" aria-hidden="true">
        <path className="edge-source" d="M150 95 C240 95 245 200 350 200"/>
        <path className="edge-predicted" d="M435 200 C525 200 535 95 625 95"/>
        <path className="edge-curated" d="M710 95 C795 95 805 200 875 200"/>
        <path className="edge-measured" d="M435 200 C530 200 535 320 625 320"/>
        <path className="edge-measured" d="M710 320 C790 320 795 200 875 200"/>
      </svg>
      <button className={`gs-node source ${active === 'source' ? 'active' : ''}`} onClick={() => setActive('source')}><Database/><span>NP SOURCE</span><strong>{sourceNodeLabel}</strong></button>
      <button className={`gs-node compound ${active === 'compound' ? 'active' : ''}`} onClick={() => setActive('compound')}><Atom/><span>COMPOUND</span><strong>{row.name}</strong></button>
      <button className={`gs-node target ${active === 'target' ? 'active' : ''}`} onClick={() => setActive('target')}><Target/><span>TARGET</span><strong>PTP1B</strong></button>
      <button className={`gs-node disease ${active === 'disease' ? 'active' : ''}`} onClick={() => setActive('disease')}><Dna/><span>DISEASE</span><strong>Type 2 diabetes</strong></button>
      <button className={`gs-node assay ${active === 'assay' ? 'active' : ''}`} onClick={() => setActive('assay')}><TestTube2/><span>ASSAY</span><strong>IC₅₀</strong></button>
      <button className={`gs-node publication ${active === 'publication' ? 'active' : ''}`} onClick={() => setActive('publication')}><GitBranch/><span>EVIDENCE</span><strong>{row.measurements ?? 'Provenance'}</strong></button>
    </div>
    <div className="gs-node-inspector">
      <div className={`gs-node-icon ${active}`}><Network/></div>
      <div><span>{nodeCopy[active].label} · {nodeCopy[active].edge}</span><h3>{node.title}</h3><p>{node.subtitle}</p></div>
      <div className="gs-evidence-note"><ShieldCheck/><span>{node.evidence}</span></div>
    </div>
  </section>
}

export function GraphScreeningWorkspace({ data }: { data: Dashboard }) {
  const [rows, setRows] = useState<ScreeningRow[]>([])
  const [selectedKey, setSelectedKey] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const [progress, setProgress] = useState('Ready')
  const [error, setError] = useState<string | null>(null)
  const [mode, setMode] = useState<'demo' | 'artifact'>('demo')
  const [query, setQuery] = useState('')
  const [domainFilter, setDomainFilter] = useState<'all' | 'in' | 'out'>('all')
  const [alertFreeOnly, setAlertFreeOnly] = useState(false)
  const [maxMolecularWeight, setMaxMolecularWeight] = useState(650)
  const [maxLogP, setMaxLogP] = useState(6)
  const [sortBy, setSortBy] = useState<'score' | 'potency' | 'similarity' | 'weight'>('score')
  const [pageSize, setPageSize] = useState(10)
  const [page, setPage] = useState(1)
  const [artifactMetadata, setArtifactMetadata] = useState<ScreeningArtifactResponse['metadata']>(null)
  const controller = useRef<AbortController | null>(null)

  useEffect(() => () => controller.current?.abort(), [])

  const demoCandidates = useMemo(() => {
    const verified = data.candidates.filter(candidate => candidate.smiles && isVerifiedNaturalProduct(candidate))
    return verified.length ? verified : data.candidates.filter(candidate => candidate.smiles).slice(0, 2)
  }, [data.candidates])

  const filteredRows = useMemo(() => rows.filter(row => {
    const matches = `${row.name} ${row.id} ${row.source}`.toLowerCase().includes(query.trim().toLowerCase())
    return matches
      && (domainFilter === 'all' || (domainFilter === 'in' ? row.prediction.applicability_domain.inside : !row.prediction.applicability_domain.inside))
      && (!alertFreeOnly || row.alerts.length === 0)
      && row.descriptors.molecular_weight <= maxMolecularWeight
      && row.descriptors.logp <= maxLogP
  }).sort((left, right) => {
    if (sortBy === 'potency') return right.prediction.predicted_pchembl - left.prediction.predicted_pchembl
    if (sortBy === 'similarity') return right.prediction.applicability_domain.nearest_training_tanimoto - left.prediction.applicability_domain.nearest_training_tanimoto
    if (sortBy === 'weight') return left.descriptors.molecular_weight - right.descriptors.molecular_weight
    return right.score - left.score
  }), [alertFreeOnly, domainFilter, maxLogP, maxMolecularWeight, query, rows, sortBy])
  const pageCount = Math.max(1, Math.ceil(filteredRows.length / pageSize))
  const visibleRows = filteredRows.slice((page - 1) * pageSize, page * pageSize)
  const selected = visibleRows.find(row => row.key === selectedKey) ?? visibleRows[0] ?? null

  useEffect(() => setPage(1), [alertFreeOnly, domainFilter, maxLogP, maxMolecularWeight, pageSize, query, sortBy])
  useEffect(() => {
    if (page > pageCount) setPage(pageCount)
  }, [page, pageCount])

  async function runDemo() {
    controller.current?.abort()
    const activeController = new AbortController()
    controller.current = activeController
    setMode('demo')
    setRows([])
    setSelectedKey(null)
    setError(null)
    setArtifactMetadata(null)
    setLoading(true)
    const output: ScreeningRow[] = []
    try {
      for (let index = 0; index < demoCandidates.length; index += 1) {
        const candidate = demoCandidates[index]
        setProgress(`Predicting ${index + 1}/${demoCandidates.length} · ${candidate.name}`)
        const prediction = await predictMolecule(candidate.smiles!, activeController.signal)
        output.push({
          rank: 0,
          key: candidate.id,
          id: candidate.id,
          name: candidate.name,
          smiles: candidate.smiles!,
          source: candidate.natural_product_status,
          organisms: candidate.organism_count > 0 ? `${candidate.organism_count} reported organism context(s)` : 'Source organism not joined in this snapshot',
          evidence: `${candidate.measurements} ChEMBL measurements · ${candidate.documents} document(s)`,
          measurements: candidate.measurements,
          assays: candidate.independent_assays,
          score: demoScore(prediction),
          prediction,
          alerts: [],
          descriptors: prediction.descriptors,
        })
      }
      const ranked = output.sort((left, right) => right.score - left.score).map((row, index) => ({ ...row, rank: index + 1 }))
      setRows(ranked)
      setSelectedKey(ranked[0]?.key ?? null)
      setProgress('Complete')
    } catch (cause) {
      if (cause instanceof DOMException && cause.name === 'AbortError') return
      setError(cause instanceof Error ? cause.message : 'The live screening demonstration failed.')
    } finally {
      setLoading(false)
      controller.current = null
    }
  }

  async function loadArtifact() {
    controller.current?.abort()
    const activeController = new AbortController()
    controller.current = activeController
    setMode('artifact')
    setRows([])
    setSelectedKey(null)
    setError(null)
    setLoading(true)
    setProgress('Loading real NPASS shortlist')
    try {
      const result = await getScreeningCandidates(250, activeController.signal)
      const loaded = result.candidates.map(artifactRow).filter((row): row is ScreeningRow => row !== null)
      setRows(loaded)
      setArtifactMetadata(result.metadata)
      setSelectedKey(loaded[0]?.key ?? null)
      setProgress('Complete')
    } catch (cause) {
      if (cause instanceof DOMException && cause.name === 'AbortError') return
      if (cause instanceof ApiError && cause.status === 503) {
        setError('No real screening artifact is configured. Run the demo now, or generate ranked_candidates.jsonl and restart the API with PHYTOTRUTH_SCREENING_RESULTS.')
      } else setError(cause instanceof Error ? cause.message : 'The screening artifact could not be loaded.')
    } finally {
      setLoading(false)
      controller.current = null
    }
  }

  const inDomainCount = rows.filter(row => row.prediction.applicability_domain.inside).length
  const bestScore = rows.length ? Math.max(...rows.map(row => row.score)).toFixed(2) : '—'
  const screenedCount = mode === 'artifact' && artifactMetadata ? artifactMetadata.screened_compounds : rows.length

  return <>
    <header className="page-title gs-title">
      <div><span className="section-kicker">KNOWLEDGE GRAPH · MULTI-TARGET SCREENING</span><h1>From natural-product provenance to an auditable shortlist</h1><p>Screen verified structures, inspect uncertainty, and follow every candidate through source, target, disease, assay and evidence relationships.</p></div>
      <div className="gs-live-chip"><span/><div><strong>PTP1B · CHEMBL335</strong><small>Type 2 diabetes research context</small></div></div>
    </header>

    <section className="panel gs-control-panel">
      <div className="gs-control-copy"><div className="gs-control-icon"><Sparkles/></div><div><span className="section-kicker">TEST SCENARIO</span><h2>Prioritize natural products for PTP1B follow-up</h2><p>The live demo runs real RDKit + XGBoost inference on the verified natural products in the bundled evidence slice. Real shortlist mode reads a precomputed, checksum-tracked screening artifact and reports its exact scope.</p></div></div>
      <div className="gs-selectors">
        <label><span>Disease</span><select aria-label="Screening disease" defaultValue="MONDO:0005148"><option value="MONDO:0005148">Type 2 diabetes mellitus</option></select></label>
        <label><span>Target</span><select aria-label="Screening target" defaultValue="CHEMBL335"><option value="CHEMBL335">PTP1B · CHEMBL335</option></select></label>
        <label><span>Source scope</span><select aria-label="Natural product source" defaultValue="verified"><option value="verified">Verified natural-product structures</option></select></label>
      </div>
      <div className="gs-actions">
        <button className="gs-run" onClick={() => void runDemo()} disabled={loading}>{loading && mode === 'demo' ? <LoaderCircle className="spin"/> : <Play/>} Run live demo <span>{demoCandidates.length} compounds</span></button>
        <button className="gs-load" onClick={() => void loadArtifact()} disabled={loading}><Database/> Load real shortlist</button>
        <div className="gs-action-note"><ShieldCheck/><span>Measured evidence and ML predictions remain visibly separate.</span></div>
      </div>
      {loading && <div className="gs-progress" role="status"><i/><span>{progress}</span></div>}
      {error && <div className="prediction-error gs-error" role="alert"><CircleAlert/><span>{error}</span>{mode === 'artifact' && <button onClick={() => void runDemo()}>Run demo instead</button>}</div>}
    </section>

    {rows.length === 0 && !loading && <section className="gs-empty panel">
      <div className="gs-empty-visual"><Network/><i/><i/><i/></div>
      <div><span className="section-kicker">READY TO SCREEN</span><h2>Run the example to create a live shortlist</h2><p>You will see real model predictions, applicability-domain checks, transparent scoring and an interactive semantic graph.</p><button onClick={() => void runDemo()}><Play/> Start guided example</button></div>
      <ol><li><b>01</b><span>Predict target activity</span></li><li><b>02</b><span>Penalize uncertainty</span></li><li><b>03</b><span>Rank candidates</span></li><li><b>04</b><span>Audit graph evidence</span></li></ol>
    </section>}

    {rows.length > 0 && <>
      <div className="gs-stats">
        <Stat label="Screened" value={screenedCount.toLocaleString()} note={mode === 'demo' ? 'live traceable slice' : `${artifactMetadata?.source_filter || 'NP'} real structures`}/>
        <Stat label="In domain" value={`${inDomainCount}/${rows.length}`} note="meets similarity threshold"/>
        <Stat label="Best score" value={bestScore} note="uncertainty-adjusted"/>
        <Stat label="Graph path" value="6 nodes" note="source to disease evidence"/>
      </div>

      <div className="gs-results-layout">
        <section className="panel gs-results-panel">
          <div className="panel-head gs-panel-head"><div><span className="section-kicker">RANKED CANDIDATES</span><h2>Experimental follow-up queue</h2></div><span className="gs-mode-pill">{mode === 'demo' ? 'LIVE DEMO' : 'FULL ARTIFACT'}</span></div>
          {mode === 'artifact' && artifactMetadata && <div className="gs-artifact-banner"><Database/><div><strong>Real NPASS screening subset</strong><span>{artifactMetadata.screened_compounds.toLocaleString()} structures evaluated · {artifactMetadata.property_filtered_compounds.toLocaleString()} property-filtered · {artifactMetadata.in_domain_candidates} IN + {artifactMetadata.out_of_domain_candidates} OUT retained</span></div><code>{artifactMetadata.natural_product_database_sha256.slice(0, 10)}…</code></div>}
          <div className="gs-filterbar">
            <div className="search"><Search/><input aria-label="Search screening candidates" value={query} onChange={event => setQuery(event.target.value)} placeholder="Search compound or source"/></div>
            <label className="gs-domain-filter"><span>Domain</span><select aria-label="Applicability domain" value={domainFilter} onChange={event => setDomainFilter(event.target.value as typeof domainFilter)}><option value="all">All · {rows.length}</option><option value="in">IN · {inDomainCount}</option><option value="out">OUT · {rows.length - inDomainCount}</option></select></label>
            <label><input type="checkbox" checked={alertFreeOnly} onChange={event => setAlertFreeOnly(event.target.checked)}/><span>No alerts</span></label>
          </div>
          <div className="gs-advanced-filters">
            <label><span>Sort</span><select aria-label="Sort shortlist" value={sortBy} onChange={event => setSortBy(event.target.value as typeof sortBy)}><option value="score">Priority score</option><option value="potency">Predicted pChEMBL</option><option value="similarity">Tanimoto similarity</option><option value="weight">Molecular weight</option></select></label>
            <label><span>Max MW</span><select aria-label="Maximum molecular weight" value={maxMolecularWeight} onChange={event => setMaxMolecularWeight(Number(event.target.value))}><option value="350">350</option><option value="500">500</option><option value="650">650</option></select></label>
            <label><span>Max LogP</span><select aria-label="Maximum LogP" value={maxLogP} onChange={event => setMaxLogP(Number(event.target.value))}><option value="3">3</option><option value="5">5</option><option value="6">6</option></select></label>
            <label><span>Rows</span><select aria-label="Rows per page" value={pageSize} onChange={event => setPageSize(Number(event.target.value))}><option value="10">10</option><option value="25">25</option><option value="50">50</option><option value="100">100</option></select></label>
          </div>
          <div className="gs-candidate-list">
            {visibleRows.map(row => <button key={row.key} className={`gs-candidate ${selected?.key === row.key ? 'active' : ''}`} onClick={() => setSelectedKey(row.key)}>
              <span className="gs-rank">{row.rank.toString().padStart(2, '0')}</span>
              <span className="gs-candidate-atom"><Atom/></span>
              <span className="gs-candidate-name"><strong>{row.name}</strong><small>{row.id}</small></span>
              <span className="gs-score"><small>priority</small><strong>{row.score.toFixed(2)}</strong></span>
              <span className={`gs-domain ${row.prediction.applicability_domain.inside ? 'inside' : 'outside'}`}>{row.prediction.applicability_domain.inside ? <CheckCircle2/> : <AlertTriangle/>}{row.prediction.applicability_domain.inside ? 'IN' : 'OUT'}</span>
              <ChevronRight/>
            </button>)}
            {filteredRows.length === 0 && <div className="gs-no-results"><Filter/><strong>No candidates pass these filters</strong><button onClick={() => { setQuery(''); setDomainFilter('all'); setAlertFreeOnly(false) }}>Reset filters</button></div>}
          </div>
          <div className="panel-foot gs-pagination"><span><FlaskConical/> {filteredRows.length} candidates match</span><div><button onClick={() => setPage(value => Math.max(1, value - 1))} disabled={page === 1}>Previous</button><span>Page {page} / {pageCount}</span><button onClick={() => setPage(value => Math.min(pageCount, value + 1))} disabled={page === pageCount}>Next</button></div></div>
        </section>

        {selected && <aside className="panel gs-detail-panel">
          <div className="gs-detail-head"><div><span className="section-kicker">SELECTED CANDIDATE</span><h2>{selected.name}</h2><p>{selected.source}</p></div><span className={`result-domain-pill ${selected.prediction.applicability_domain.inside ? 'inside' : 'outside'}`}>{selected.prediction.applicability_domain.inside ? <CheckCircle2/> : <AlertTriangle/>}{selected.prediction.applicability_domain.inside ? 'In domain' : 'Review required'}</span></div>
          <div className="gs-molecule"><MoleculeViewer smiles={selected.smiles}/></div>
          <div className="gs-detail-values"><div><span>Predicted pChEMBL</span><strong>{selected.prediction.predicted_pchembl.toFixed(3)}</strong><small>ML prediction</small></div><div><span>Priority score</span><strong>{selected.score.toFixed(3)}</strong><small>uncertainty-adjusted</small></div></div>
          <div className="gs-interval"><span>90% residual interval</span><strong>{selected.prediction.prediction_interval.lower.toFixed(2)} — {selected.prediction.prediction_interval.upper.toFixed(2)}</strong></div>
          <div className="gs-similarity"><div><span>Tanimoto similarity</span><strong>{selected.prediction.applicability_domain.nearest_training_tanimoto.toFixed(3)}</strong></div><div className="gs-similarity-bar"><i style={{ width: `${selected.prediction.applicability_domain.nearest_training_tanimoto * 100}%` }}/><b style={{ left: `${selected.prediction.applicability_domain.threshold * 100}%` }}/></div><small>Domain threshold {selected.prediction.applicability_domain.threshold.toFixed(3)}</small></div>
          {!selected.prediction.applicability_domain.inside && <div className="gs-alert-box"><AlertTriangle/><div><strong>Out-of-domain prediction — review required</strong><p>This chemistry is below the training-similarity threshold. Treat its predicted activity as exploratory and require stronger experimental validation.</p></div></div>}
          <div className="gs-alert-box"><CircleAlert/><div><strong>{selected.alerts.length ? `${selected.alerts.length} structural alert(s)` : 'No stored structural alerts'}</strong><p>{selected.alerts.length ? selected.alerts.join(' · ') : mode === 'demo' ? 'The live demo shows model/domain behavior; full PAINS/Brenk scanning is recorded by the batch artifact.' : 'PAINS/Brenk scan completed in the batch pipeline.'}</p></div></div>
        </aside>}
      </div>
      {selected && <EvidenceGraph row={selected}/>} 
      <div className="gs-integrity"><ShieldCheck/><div><strong>Research integrity boundary</strong><p>Ranking means “prioritize for experimental validation.” It does not establish efficacy, safety, toxicity, dosing, causality or clinical benefit.</p></div><a href={apiDocsUrl()} target="_blank" rel="noreferrer">Open API contract <ExternalLink/></a></div>
    </>}
  </>
}
