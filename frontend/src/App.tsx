import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  Activity, ArrowDown, ArrowRight, ArrowUp, Atom, BookOpen, Check,
  ChevronRight, CircleAlert, Database, Dna, ExternalLink, FlaskConical,
  BrainCircuit, GitCompareArrows, Leaf, Menu, Microscope, Network, Search, ShieldCheck, Sparkles, X,
} from 'lucide-react'
import {
  CartesianGrid, Cell, ResponsiveContainer, Scatter, ScatterChart, Tooltip,
  XAxis, YAxis,
} from 'recharts'
import { fallbackDashboard } from './demo'
import { getDashboard, moleculeSvgUrl } from './api'
import { ExperimentWorkspace } from './ExperimentWorkspace'
import { GraphScreeningWorkspace } from './GraphScreeningWorkspace'
import { PredictionWorkspace } from './PredictionWorkspace'
import { ScientificValidationWorkspace } from './ScientificValidationWorkspace'
import { StudyPlannerWorkspace } from './StudyPlannerWorkspace'
import type { Candidate, Dashboard } from './types'

type View = 'Discovery' | 'Graph & Screening' | 'Predict' | 'Evidence' | 'Scientific Validation' | 'Study Planner' | 'Reversal Lab' | 'Data Sources'

const nav: { label: View; icon: typeof Atom }[] = [
  { label: 'Discovery', icon: Atom },
  { label: 'Graph & Screening', icon: Network },
  { label: 'Predict', icon: Sparkles },
  { label: 'Evidence', icon: ShieldCheck },
  { label: 'Study Planner', icon: BrainCircuit },
  { label: 'Scientific Validation', icon: Microscope },
  { label: 'Reversal Lab', icon: GitCompareArrows },
  { label: 'Data Sources', icon: Database },
]

function scoreTone(value: number) {
  return value >= .65 ? 'mint' : value >= .4 ? 'amber' : 'rose'
}

function MoleculeArt({ candidate }: { candidate: Candidate }) {
  const [failed, setFailed] = useState(false)
  if (!candidate.smiles || failed) {
    return <div className="molecule-fallback"><Atom size={44}/><span>Structure preview unavailable</span></div>
  }
  return <img
    className="molecule-img"
    src={moleculeSvgUrl(candidate.smiles)}
    alt={`RDKit 2D structure for ${candidate.name}`}
    onError={() => setFailed(true)}
  />
}

function RankDelta({ value }: { value: number }) {
  if (!value) return <span className="delta neutral">—</span>
  return <span className={`delta ${value > 0 ? 'up' : 'down'}`}>
    {value > 0 ? <ArrowUp size={12}/> : <ArrowDown size={12}/>} {Math.abs(value)}
  </span>
}

function Metric({ label, value, sub, icon: Icon }: { label: string; value: string | number; sub: string; icon: typeof Atom }) {
  return <div className="metric-card">
    <div className="metric-icon"><Icon size={18}/></div>
    <div><span>{label}</span><strong>{value}</strong><small>{sub}</small></div>
  </div>
}

function CandidateDrawer({ candidate, close }: { candidate: Candidate; close: () => void }) {
  const closeButton = useRef<HTMLButtonElement>(null)

  useEffect(() => {
    const previousFocus = document.activeElement instanceof HTMLElement ? document.activeElement : null
    const previousOverflow = document.body.style.overflow
    document.body.style.overflow = 'hidden'
    closeButton.current?.focus()
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') close()
    }
    window.addEventListener('keydown', handleKeyDown)
    return () => {
      window.removeEventListener('keydown', handleKeyDown)
      document.body.style.overflow = previousOverflow
      previousFocus?.focus()
    }
  }, [close])

  return <div className="drawer-wrap" role="dialog" aria-modal="true" aria-label={`${candidate.name} evidence details`}>
    <button className="drawer-scrim" onClick={close} aria-label="Close details"/>
    <aside className="drawer">
      <button ref={closeButton} className="icon-btn drawer-close" onClick={close} aria-label="Close evidence details"><X size={18}/></button>
      <div className="eyebrow">Evidence dossier · {candidate.id}</div>
      <h2>{candidate.name}</h2>
      <p className="muted">{candidate.natural_product_status}</p>
      <div className="drawer-molecule"><MoleculeArt candidate={candidate}/></div>
      <div className="score-hero">
        <div className={`score-ring ${scoreTone(candidate.reliability)}`} style={{'--score': `${candidate.reliability * 360}deg`} as React.CSSProperties}>
          <div><strong>{Math.round(candidate.reliability * 100)}</strong><span>reliability</span></div>
        </div>
        <div>
          <span className="micro">Conservative evidence score</span>
          <div className="score-number">{candidate.evidence_score.toFixed(2)}</div>
          <p>pChEMBL median {candidate.median_pchembl.toFixed(2)} · 95% bootstrap interval {candidate.confidence_interval[0]}–{candidate.confidence_interval[1]}</p>
        </div>
      </div>
      <div className="detail-grid">
        <div><span>Measurements</span><strong>{candidate.measurements}</strong></div>
        <div><span>Independent assays</span><strong>{candidate.independent_assays}</strong></div>
        <div><span>Documents</span><strong>{candidate.documents}</strong></div>
        <div><span>pChEMBL range</span><strong>{candidate.pchembl_range.toFixed(2)}</strong></div>
      </div>
      <section className="drawer-section">
        <h3>Truth signals</h3>
        {candidate.flags.length ? candidate.flags.map(flag => <div className="flag-row" key={flag}><CircleAlert size={15}/><span>{flag}</span></div>) : <div className="flag-row clean"><Check size={15}/><span>No automated inconsistency flags in this slice</span></div>}
      </section>
      <section className="drawer-section">
        <h3>How to read this</h3>
        <p>The reliability number describes agreement and evidence support in this dataset slice. It is not a probability of efficacy, safety, or clinical success.</p>
      </section>
      <a className="source-link" href={`https://www.ebi.ac.uk/chembl/explore/compound/${candidate.id}`} target="_blank" rel="noreferrer">
        Open original ChEMBL record <ExternalLink size={14}/>
      </a>
    </aside>
  </div>
}

function Discovery({ data, select }: { data: Dashboard; select: (c: Candidate) => void }) {
  const [query, setQuery] = useState('')
  const normalizedQuery = query.trim().toLowerCase()
  const candidates = data.candidates.filter(item => `${item.name} ${item.id}`.toLowerCase().includes(normalizedQuery))
  return <>
    <section className="hero">
      <div>
        <div className="eyebrow"><span className="live-dot"/> LIVE EVIDENCE WORKSPACE · ETHNOHACK 2026</div>
        <h1>Find the signal.<br/><em>Question the evidence.</em></h1>
        <p>PhytoTruth surfaces promising natural products, then stress-tests every ranking against assay disagreement, provenance gaps, and uncertainty.</p>
      </div>
      <div className="target-card">
        <span className="micro">ACTIVE TARGET</span>
        <div className="target-title"><div><Dna size={22}/></div><strong>PTP1B</strong><span>CHEMBL335</span></div>
        <p>{data.metadata.target_name}</p>
        <div className="target-footer"><span>Homo sapiens</span><span>IC₅₀ evidence</span></div>
      </div>
    </section>

    <div className="metrics">
      <Metric label="Candidate slice" value={data.summary.molecules} sub="traceable molecules" icon={Leaf}/>
      <Metric label="Evidence points" value={data.summary.measurements} sub="measured activities" icon={Activity}/>
      <Metric label="Independent assays" value={data.summary.independent_assays} sub="distinct protocols" icon={FlaskConical}/>
      <Metric label="Truth alerts" value={data.summary.flagged_candidates} sub="need researcher review" icon={CircleAlert}/>
    </div>

    <section className="panel candidates-panel">
      <div className="panel-head">
        <div><span className="section-kicker">EVIDENCE-AWARE RANKING</span><h2>Candidate frontier</h2></div>
        <div className="search"><Search size={15}/><input aria-label="Search candidates" value={query} onChange={event => setQuery(event.target.value)} placeholder="Search compound or ChEMBL ID"/></div>
      </div>
      <div className="table-scroll">
        <table>
          <thead><tr><th>#</th><th>Candidate</th><th>Evidence score</th><th>Reliability</th><th>Assays</th><th>Spread</th><th>vs potency</th><th></th></tr></thead>
          <tbody>{candidates.map(candidate => <tr key={candidate.id} onClick={() => select(candidate)}>
            <td><span className="rank">{candidate.evidence_rank.toString().padStart(2,'0')}</span></td>
            <td><div className="candidate-name"><div className="mini-atom"><Atom size={17}/></div><div><strong>{candidate.name}</strong><span>{candidate.id}</span></div></div></td>
            <td><strong className="tabular">{candidate.evidence_score.toFixed(2)}</strong><span className="tiny-label">conservative pChEMBL</span></td>
            <td><div className="reliability"><div className="bar"><i className={scoreTone(candidate.reliability)} style={{width:`${candidate.reliability*100}%`}}/></div><span>{Math.round(candidate.reliability*100)}%</span></div></td>
            <td className="tabular">{candidate.independent_assays}</td>
            <td><span className={candidate.pchembl_range >= 1.5 ? 'warning-pill' : 'plain-pill'}>{candidate.pchembl_range.toFixed(2)}</span></td>
            <td><RankDelta value={candidate.rank_change}/></td>
            <td><button className="row-open" aria-label={`Open ${candidate.name}`}><ChevronRight size={17}/></button></td>
          </tr>)}{candidates.length === 0 && <tr className="empty-row"><td colSpan={8}><Search size={20}/><strong>No candidates match “{query.trim()}”</strong><span>Try a compound name or ChEMBL identifier.</span><button type="button" onClick={() => setQuery('')}>Clear search</button></td></tr>}</tbody>
        </table>
      </div>
      <div className="panel-foot"><span><ShieldCheck size={14}/> Scores penalize inconsistency and sparse support</span><span>Source snapshot: {data.metadata.retrieved_at}</span></div>
    </section>
  </>
}

function EvidenceView({ data, select }: { data: Dashboard; select: (c: Candidate) => void }) {
  const chart = data.candidates.map(item => ({...item, potency:item.max_pchembl, evidence:item.evidence_score}))
  return <>
    <header className="page-title"><div><span className="section-kicker">EVIDENCE MAP</span><h1>Where potency and truth diverge</h1><p>Each point is a molecule. Distance below the diagonal is the penalty for inconsistent or weakly supported evidence.</p></div></header>
    <div className="evidence-layout">
      <section className="panel chart-panel">
        <div className="panel-head"><div><h2>Potency claim vs evidence score</h2><p className="muted">Click a point to inspect the audit trail</p></div><span className="legend"><i/> reliability</span></div>
        <div className="chart-wrap"><ResponsiveContainer width="100%" height="100%"><ScatterChart margin={{top:18,right:18,bottom:12,left:0}}>
          <CartesianGrid strokeDasharray="3 5" stroke="#1b3430"/>
          <XAxis dataKey="potency" name="max potency" type="number" domain={[4,10]} tick={{fill:'#7e9992',fontSize: 12}} axisLine={{stroke:'#29443e'}} label={{value:'Maximum reported pChEMBL',position:'insideBottom',offset:-4,fill:'#7e9992',fontSize: 12}}/>
          <YAxis dataKey="evidence" name="evidence score" type="number" domain={[4,10]} tick={{fill:'#7e9992',fontSize: 12}} axisLine={{stroke:'#29443e'}}/>
          <Tooltip cursor={{strokeDasharray:'3 3'}} content={({active,payload}) => active && payload?.[0] ? <div className="chart-tip"><strong>{payload[0].payload.name}</strong><span>Max {payload[0].payload.potency} → evidence {payload[0].payload.evidence}</span></div> : null}/>
          <Scatter data={chart} onClick={(point) => select(point as unknown as Candidate)}>{chart.map((item) => <Cell key={item.id} fill={item.reliability >= .65 ? '#6bf5bb' : item.reliability >= .4 ? '#f0ba69' : '#ff7185'} stroke="#07110f" strokeWidth={3}/>)}</Scatter>
        </ScatterChart></ResponsiveContainer></div>
      </section>
      <aside className="panel truth-panel">
        <span className="section-kicker">AUTOMATED REVIEW QUEUE</span><h2>Signals needing attention</h2>
        {data.candidates.filter(c => c.flags.length).map(candidate => <button className="truth-item" key={candidate.id} onClick={() => select(candidate)}>
          <div className={`truth-icon ${scoreTone(candidate.reliability)}`}><CircleAlert size={16}/></div>
          <div><strong>{candidate.name}</strong><span>{candidate.flags.join(' · ')}</span></div><ChevronRight size={15}/>
        </button>)}
        <div className="method-note"><ShieldCheck size={17}/><div><strong>Transparent by design</strong><p>Every flag comes from a deterministic rule in <code>evidence.py</code>, not an opaque language-model judgment.</p></div></div>
      </aside>
    </div>
  </>
}

function ReversalLab({ data }: { data: Dashboard }) {
  const potency = [...data.candidates].sort((a,b)=>a.potency_rank-b.potency_rank)
  const evidence = [...data.candidates].sort((a,b)=>a.evidence_rank-b.evidence_rank)
  return <>
    <header className="page-title"><div><span className="section-kicker">EVIDENCE REVERSAL EXPERIMENT</span><h1>Would yesterday’s ranking survive tomorrow’s assay?</h1><p>Compare a best-result potency ranking with a reliability-aware ranking. Production experiments use compound-held-out evidence and scaffold-grouped model validation.</p></div><div className="lab-status"><span/><div><strong>Protocol ready</strong><small>requires full dataset import</small></div></div></header>
    <div className="experiment-banner"><Sparkles size={20}/><div><strong>The hypothesis</strong><p>Penalizing contradictions should improve recovery of compounds whose activity replicates in held-out assays—even when their single best result is less spectacular.</p></div></div>
    <div className="rank-compare">
      <section className="panel rank-column"><div className="rank-head"><div><span className="bad-dot"/>BASELINE</div><h2>Potency-only</h2><p>Rank by maximum reported pChEMBL</p></div>{potency.map(c=><div className="rank-item" key={c.id}><b>{c.potency_rank}</b><div><strong>{c.name}</strong><span>best pChEMBL {c.max_pchembl.toFixed(2)}</span></div>{c.flags.length > 0 && <CircleAlert className="amber-text" size={16}/>}</div>)}</section>
      <div className="versus"><ArrowRight size={18}/><span>RE-RANK</span></div>
      <section className="panel rank-column preferred"><div className="rank-head"><div><span className="good-dot"/>PHYTOTRUTH</div><h2>Evidence-aware</h2><p>Central potency − inconsistency + support</p></div>{evidence.map(c=><div className="rank-item" key={c.id}><b>{c.evidence_rank}</b><div><strong>{c.name}</strong><span>evidence {c.evidence_score.toFixed(2)} · reliability {Math.round(c.reliability*100)}%</span></div><RankDelta value={c.rank_change}/></div>)}</section>
    </div>
    <section className="panel protocol"><div><span className="section-kicker">LEAKAGE BARRIERS</span><h2>Experiment contract</h2></div><div className="protocol-steps">
      <div><span>01</span><strong>Group by molecule</strong><p>Replicates never cross train and held-out evidence.</p></div>
      <div><span>02</span><strong>Split by scaffold</strong><p>Murcko groups prevent close analog leakage.</p></div>
      <div><span>03</span><strong>Freeze ranking</strong><p>Only training evidence determines priority.</p></div>
      <div><span>04</span><strong>Score held-out</strong><p>nDCG@k, EF@k, precision@k and Kendall τ.</p></div>
    </div></section>
  </>
}

function DataSources({ data }: { data: Dashboard }) {
  const sources = [
    {name:'ChEMBL', status:'Connected snapshot', count:`${data.summary.measurements} activity records`, text:'Measured bioactivities, assay and document identifiers.', color:'green', url:'https://www.ebi.ac.uk/chembl/'},
    {name:'COCONUT 2.0', status:'Parser verified', count:'1,000-row real-file check', text:'Natural-product structures from the official 2024-08 Drug Discovery export; full-release ingestion remains pending.', color:'cyan', url:'https://coconut.naturalproducts.net/download'},
    {name:'NPASS', status:'Release indexed', count:'198,028 unique structures', text:'Complete NPASS 3.0 structure export indexed; activity and organism files are not yet joined.', color:'violet', url:'https://bidd.group/NPASS/'},
  ]
  return <>
    <header className="page-title"><div><span className="section-kicker">DATA LINEAGE</span><h1>Every conclusion starts with a source</h1><p>Raw identifiers are retained end-to-end. Imports are versioned; predictions never overwrite experiments.</p></div></header>
    <div className="source-grid">{sources.map(source=><a className="source-card" href={source.url} target="_blank" rel="noreferrer" key={source.name}>
      <div className={`source-mark ${source.color}`}><Database size={22}/></div><div className="source-status"><i/> {source.status}</div><h2>{source.name}</h2><strong>{source.count}</strong><p>{source.text}</p><span>Open official source <ExternalLink size={13}/></span>
    </a>)}</div>
    <section className="panel lineage"><div><span className="section-kicker">PROVENANCE PIPELINE</span><h2>No silent transformations</h2></div><div className="lineage-flow"><div><Database/><strong>Raw source</strong><span>ID + release + row</span></div><ArrowRight/><div><FlaskConical/><strong>Harmonize</strong><span>units + relations</span></div><ArrowRight/><div><ShieldCheck/><strong>Audit evidence</strong><span>flags + uncertainty</span></div><ArrowRight/><div><Atom/><strong>Rank</strong><span>measured ≠ predicted</span></div></div></section>
    <div className="integrity-note"><BookOpen size={20}/><div><strong>Scientific integrity boundary</strong><p>The bundled dashboard is a small real ChEMBL slice for demonstration. It does not claim the unlabeled molecules are natural products. Run COCONUT/NPASS ingestion and structure cross-matching before a natural-product screen.</p></div></div>
  </>
}

export default function App() {
  const [view, setView] = useState<View>('Discovery')
  const [data, setData] = useState<Dashboard>(fallbackDashboard)
  const [source, setSource] = useState<'loading'|'api'|'snapshot'>('loading')
  const [loadError, setLoadError] = useState<string | null>(null)
  const [selected, setSelected] = useState<Candidate | null>(null)
  const [mobileNav, setMobileNav] = useState(false)

  const loadDashboard = useCallback((signal?: AbortSignal) => {
    setSource('loading')
    setLoadError(null)
    return getDashboard(signal)
      .then(value => { setData(value); setSource('api') })
      .catch(error => {
        if (error instanceof DOMException && error.name === 'AbortError') return
        setData(fallbackDashboard)
        setSource('snapshot')
        setLoadError(error instanceof Error ? error.message : 'The live API could not be reached.')
      })
  }, [])

  useEffect(() => {
    const controller = new AbortController()
    void loadDashboard(controller.signal)
    return () => controller.abort()
  }, [loadDashboard])

  const content = useMemo(() => {
    if (view === 'Graph & Screening') return <GraphScreeningWorkspace data={data}/>
    if (view === 'Predict') return <PredictionWorkspace candidates={data.candidates}/>
    if (view === 'Evidence') return <EvidenceView data={data} select={setSelected}/>
    if (view === 'Scientific Validation') return <ScientificValidationWorkspace/>
    if (view === 'Study Planner') return <StudyPlannerWorkspace/>
    if (view === 'Reversal Lab') return <ExperimentWorkspace/>
    if (view === 'Data Sources') return <DataSources data={data}/>
    return <Discovery data={data} select={setSelected}/>
  }, [data, view])

  return <div className="app-shell">
    <aside className={`sidebar ${mobileNav ? 'open' : ''}`}>
      <div className="brand"><div className="brand-mark"><Leaf/><span/></div><div><strong>PhytoTruth</strong><span>AI</span></div></div>
      <button className="nav-close" onClick={()=>setMobileNav(false)}><X/></button>
      <nav>{nav.map(item => <button key={item.label} className={view===item.label?'active':''} onClick={()=>{setView(item.label);setMobileNav(false)}}><item.icon size={18}/><span>{item.label}</span>{item.label==='Evidence' && <b>{data.summary.flagged_candidates}</b>}</button>)}</nav>
      <div className="sidebar-bottom">
        <div className="mode-card"><Microscope size={17}/><div><strong>Research mode</strong><span>No clinical claims</span></div></div>
        <div className="version">v0.2 research MVP <span>·</span> EthnoHack ’26</div>
      </div>
    </aside>
    {mobileNav && <button className="mobile-scrim" onClick={()=>setMobileNav(false)}/>}
    <main>
      <div className="topbar"><button className="mobile-menu" onClick={()=>setMobileNav(true)} aria-label="Open navigation"><Menu/></button><div className="breadcrumb"><span>WORKSPACE</span><ChevronRight size={13}/><strong>{view}</strong></div><div className="api-state" role="status"><i className={source}/>{source === 'api' ? 'Live API' : source === 'loading' ? 'Connecting…' : 'Verified snapshot'}</div></div>
      <div className="content">
        {source === 'snapshot' && <div className="api-notice" role="alert"><CircleAlert size={17}/><div><strong>Live API unavailable — showing the bundled snapshot</strong><span>{loadError} Predictions and live structure rendering require the API.</span></div><button type="button" onClick={() => void loadDashboard()}>Retry connection</button></div>}
        {content}
      </div>
      <footer><span>For research prioritization only — not medical advice.</span><span>Source IDs and uncertainty remain visible by design.</span></footer>
    </main>
    {selected && <CandidateDrawer candidate={selected} close={()=>setSelected(null)}/>} 
  </div>
}

