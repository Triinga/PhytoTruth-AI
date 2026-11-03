import { FormEvent, useEffect, useRef, useState } from 'react'
import {
  AlertTriangle, Atom, CheckCircle2, ChevronRight, CircleAlert, FlaskConical,
  Gauge, LoaderCircle, RotateCcw, Send, ShieldAlert, Sparkles,
} from 'lucide-react'
import { ApiError, predictMolecule } from './api'
import { MoleculeViewer } from './MoleculeViewer'
import type { Candidate, MolecularDescriptors, PredictionResponse } from './types'

const descriptorLabels: Record<keyof MolecularDescriptors, string> = {
  molecular_weight: 'Molecular weight',
  logp: 'LogP',
  tpsa: 'TPSA',
  hbd: 'H-bond donors',
  hba: 'H-bond acceptors',
  rotatable_bonds: 'Rotatable bonds',
  fraction_csp3: 'Fraction Csp³',
  ring_count: 'Ring count',
  heavy_atoms: 'Heavy atoms',
  formal_charge: 'Formal charge',
}

type Example = { name: string; note: string; smiles: string }

function examplesFrom(candidates: Candidate[]): Example[] {
  const examples = candidates.filter(candidate => candidate.smiles).slice(0, 2).map(candidate => ({
    name: candidate.name,
    note: candidate.natural_product_status.includes('known') ? 'Known natural product in evidence slice' : 'Traceable ChEMBL structure',
    smiles: candidate.smiles!,
  }))
  return [...examples, { name: 'Ethanol', note: 'Small out-of-domain control', smiles: 'CCO' }]
}

export function PredictionStructure({ result }: { result: PredictionResponse }) {
  return <MoleculeViewer smiles={result.canonical_smiles}/>
}

function DomainPanel({ result }: { result: PredictionResponse }) {
  const domain = result.applicability_domain
  const ratio = Math.min(100, (domain.nearest_training_tanimoto / Math.max(domain.threshold, .001)) * 72)
  return <section className={`domain-card ${domain.inside ? 'in-domain' : 'out-domain'}`} aria-label="Applicability domain assessment">
    <div className="domain-card-head">
      {domain.inside ? <CheckCircle2/> : <ShieldAlert/>}
      <div>
        <span>APPLICABILITY DOMAIN</span>
        <strong>{domain.inside ? 'Within modeled chemistry' : 'Outside modeled chemistry'}</strong>
      </div>
    </div>
    <div className="domain-meter" aria-label={`Nearest training similarity ${domain.nearest_training_tanimoto.toFixed(3)}, threshold ${domain.threshold.toFixed(3)}`}>
      <i style={{ width: `${ratio}%` }}/><b style={{ left: '72%' }}/>
    </div>
    <div className="domain-scale"><span>Similarity {domain.nearest_training_tanimoto.toFixed(3)}</span><span>Threshold {domain.threshold.toFixed(3)}</span></div>
    <p>{domain.inside
      ? 'The structure meets the prespecified training-similarity threshold. This remains a research estimate.'
      : 'Do not treat this prediction as reliable. The molecule is less similar to training chemistry than the prespecified threshold.'}</p>
  </section>
}

export function PredictionWorkspace({ candidates }: { candidates: Candidate[] }) {
  const [smiles, setSmiles] = useState('')
  const [result, setResult] = useState<PredictionResponse | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const controller = useRef<AbortController | null>(null)
  const requestVersion = useRef(0)
  const examples = examplesFrom(candidates)

  useEffect(() => () => {
    requestVersion.current += 1
    controller.current?.abort()
  }, [])

  function invalidatePrediction() {
    requestVersion.current += 1
    controller.current?.abort()
    controller.current = null
    setLoading(false)
    setResult(null)
  }

  async function submit(event: FormEvent) {
    event.preventDefault()
    const value = smiles.trim()
    if (!value) {
      setError('Enter a SMILES structure before requesting a prediction.')
      return
    }
    if (value.length > 10_000) {
      setError('The SMILES exceeds the API limit of 10,000 characters.')
      return
    }
    controller.current?.abort()
    const activeController = new AbortController()
    controller.current = activeController
    const activeVersion = requestVersion.current + 1
    requestVersion.current = activeVersion
    setLoading(true)
    setError(null)
    setResult(null)
    try {
      const prediction = await predictMolecule(value, activeController.signal)
      if (activeVersion !== requestVersion.current) return
      setResult(prediction)
    } catch (cause) {
      if (activeVersion !== requestVersion.current) return
      if (cause instanceof DOMException && cause.name === 'AbortError') return
      if (cause instanceof ApiError) {
        if (cause.code === 'invalid_smiles' || cause.status === 422) setError('RDKit could not parse this SMILES. Check ring closures, brackets, bond symbols, and stereochemistry.')
        else if (cause.code === 'model_unavailable' || cause.status === 503) setError('The trained model is unavailable or incompatible. Check the backend model artifacts.')
        else setError(cause.message)
      } else setError('Prediction failed unexpectedly. No result was generated.')
    } finally {
      if (activeVersion === requestVersion.current) {
        setLoading(false)
        controller.current = null
      }
    }
  }

  function useExample(example: Example) {
    invalidatePrediction()
    setSmiles(example.smiles)
    setError(null)
  }

  return <>
    <header className="page-title prediction-title">
      <div><span className="section-kicker">MOLECULAR PREDICTION · REAL XGBOOST INFERENCE</span><h1>Test a structure against the model’s known chemistry</h1><p>Submit a SMILES string for a target-specific CHEMBL335 IC₅₀ estimate. The interface keeps model output, measured evidence, and experiment findings explicitly separate.</p></div>
      <div className="model-chip"><Sparkles size={17}/><div><strong>Saved model</strong><span>2,444 structures · 1,086 scaffolds</span></div></div>
    </header>

    <div className="prediction-layout">
      <section className="panel predict-form-panel">
        <div className="predict-form-head"><div className="prediction-kind"><Atom size={17}/> ML PREDICTION</div><span>Research use only</span></div>
        <form onSubmit={submit} noValidate>
          <label htmlFor="smiles-input">Molecular SMILES</label>
          <div className={`smiles-box ${error ? 'has-error' : ''}`}>
            <textarea
              id="smiles-input"
              value={smiles}
              onChange={event => { invalidatePrediction(); setSmiles(event.target.value); setError(null) }}
              placeholder="Example: CCO"
              rows={6}
              spellCheck={false}
              aria-describedby={error ? 'prediction-error' : 'smiles-help'}
              aria-invalid={Boolean(error)}
            />
            {smiles && <button type="button" className="clear-smiles" onClick={() => { invalidatePrediction(); setSmiles(''); setError(null) }} aria-label="Clear SMILES"><RotateCcw size={14}/></button>}
          </div>
          <div id="smiles-help" className="field-help"><span>{smiles.length.toLocaleString()} / 10,000 characters</span><span>Canonicalized server-side with RDKit</span></div>
          {error && <div className="prediction-error" id="prediction-error" role="alert"><CircleAlert size={16}/><span>{error}</span></div>}
          <button className="predict-button" type="submit" disabled={loading}>
            {loading ? <><LoaderCircle className="spin"/> Running RDKit + XGBoost…</> : <><Send size={16}/> Run prediction <ChevronRight size={15}/></>}
          </button>
        </form>
        <div className="example-list"><span>TRACEABLE EXAMPLES</span>{examples.map(example => <button key={example.name} type="button" onClick={() => useExample(example)}>
          <div><strong>{example.name}</strong><small>{example.note}</small></div><ChevronRight size={14}/>
        </button>)}</div>
        <div className="prediction-boundary"><FlaskConical size={17}/><p><strong>This is not an efficacy prediction.</strong> The model estimates a compatible target-specific pChEMBL context. It does not establish experimental activity, dosing, toxicity, drug efficacy, or clinical safety.</p></div>
      </section>

      <section className={`prediction-output ${result ? 'has-result' : ''}`} aria-live="polite">
        {loading && <div className="prediction-loading"><div className="orbital-loader"><Atom/><i/><i/></div><strong>Evaluating molecular features</strong><span>Canonicalizing · Morgan fingerprint · descriptors · domain check</span></div>}
        {!loading && !result && <div className="prediction-empty"><div><Gauge/></div><h2>No prediction yet</h2><p>Choose a traceable example or enter a SMILES. Results appear here with uncertainty and an applicability-domain assessment.</p></div>}
        {result && <div className="prediction-result">
          <div className="result-topline"><div><span className="prediction-kind">ML PREDICTION</span><small>{result.model.target.target_chembl_id} · {result.model.target.standard_type} · {result.model.target.target_organism}</small></div><span className={`result-domain-pill ${result.applicability_domain.inside ? 'inside' : 'outside'}`}>{result.applicability_domain.inside ? <CheckCircle2/> : <AlertTriangle/>}{result.applicability_domain.inside ? 'In domain' : 'Outside domain'}</span></div>
          <div className="result-hero">
            <PredictionStructure result={result}/>
            <div className="prediction-value"><span>Predicted pChEMBL</span><strong>{result.predicted_pchembl.toFixed(3)}</strong><small>Model estimate · not measured</small></div>
          </div>
          <div className="canonical-smiles"><span>CANONICAL SMILES</span><code>{result.canonical_smiles}</code></div>
          <div className="interval-card">
            <div><span>NOMINAL {(result.prediction_interval.confidence_level * 100).toFixed(0)}% RESIDUAL INTERVAL</span><strong>{result.prediction_interval.lower.toFixed(3)} <i>to</i> {result.prediction_interval.upper.toFixed(3)}</strong></div>
            <div className="interval-line"><i/><b style={{ left: '50%' }}/></div>
            <p>{result.prediction_interval.method}. Retrospective aggregate coverage does not guarantee molecule-specific calibration.</p>
          </div>
          <DomainPanel result={result}/>
          <section className="descriptor-section"><div className="result-section-head"><span>RDKIT DESCRIPTORS</span><small>Calculated from canonical structure</small></div><div className="descriptor-grid">
            {(Object.entries(result.descriptors) as [keyof MolecularDescriptors, number][]).map(([key, value]) => <div key={key}><span>{descriptorLabels[key]}</span><strong>{Number.isInteger(value) ? value : value.toFixed(3)}</strong></div>)}
          </div></section>
          <div className="model-provenance"><span>MODEL VERSION</span><code>{result.model.version}</code><span>FORMAT</span><code>{result.model.format}</code></div>
        </div>}
      </section>
    </div>
  </>
}
