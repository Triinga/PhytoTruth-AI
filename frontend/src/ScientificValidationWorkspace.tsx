import { useMemo, useState } from 'react'
import {
  AlertTriangle, ExternalLink, GitCommitHorizontal, Leaf, Search, ShieldAlert, X,
} from 'lucide-react'
import rawResults from './data/ethno-convergence-results.json'
import rawSensitivity from './data/ethno-convergence-sensitivity.json'
import rawPlants from './data/ethno-plants.json'
import type {
  EthnoConvergenceResults, EthnoPlant, EthnoPlantSummary, EthnoSensitivityResults,
} from './types'

const results = rawResults as unknown as EthnoConvergenceResults
const sensitivity = rawSensitivity as unknown as EthnoSensitivityResults
const summary = rawPlants as unknown as EthnoPlantSummary

const REPO = 'https://github.com/Triinga/PhytoTruth-AI/commit'

/** The public, immutable record of protocol-before-outcome. */
const timeline = [
  {
    time: '11:22:01',
    sha: '482bf0ded20eeb1d2bdc8841d9bbcae3db43b9b4',
    what: 'Protocol frozen — hypothesis, terms, statistic, seed, bootstrap, negative control, support rule',
    isOutcome: false,
  },
  {
    time: '11:25:17',
    sha: '97cb951229ca28ed495228ee7261665d188ff704',
    what: 'Implementation and the NPASS source-organism join',
    isOutcome: false,
  },
  {
    time: '13:32:39',
    sha: 'ed63264beefc5909622fb5fdd48e81e438476ef3',
    what: 'results.json — the first outcome to exist anywhere',
    isOutcome: true,
  },
]

const DEFAULT_SPECIES = 'coix lacryma-jobi'

const titleCase = (value: string) => value.replace(/\b[a-z]/g, (char) => char.toUpperCase())
const signed = (value: number, digits = 3) => `${value >= 0 ? '+' : '−'}${Math.abs(value).toFixed(digits)}`
const pct = (part: number, whole: number) => `${((100 * part) / whole).toFixed(1)}%`

/**
 * Interval plot: the primary effect with its family-clustered CI, against the
 * negative control. A point estimate is drawn as a point; only the primary
 * outcome has a prespecified bootstrap interval, so only it gets a band.
 */
function IntervalPlot() {
  const primary = results.primary_outcome
  const control = results.negative_control
  const [low, high] = primary.family_clustered_ci_95

  const domain = 0.055
  const x = (value: number) => 50 + ((value + domain) / (2 * domain)) * 600
  const ticks = [-0.04, -0.02, 0, 0.02, 0.04]

  return <figure className="interval-plot">
    <svg viewBox="0 0 700 190" role="img" aria-label={
      `Median shift in predicted pChEMBL. Traditional antidiabetic use: ${primary.median_shift_hodges_lehmann.toFixed(4)}, `
      + `95% family-clustered confidence interval ${low.toFixed(4)} to ${high.toFixed(4)}, p equals ${primary.p_value.toFixed(3)}. `
      + `Negative control, toothache use: ${control.median_shift_hodges_lehmann.toFixed(4)}, p equals ${control.p_value.toFixed(3)}. `
      + 'Both intervals include zero.'
    }>
      {ticks.map((tick) => <g key={tick}>
        <line x1={x(tick)} x2={x(tick)} y1={28} y2={132} stroke={tick === 0 ? '#5e7d74' : '#16302a'} strokeWidth={tick === 0 ? 1.5 : 1}/>
        <text x={x(tick)} y={152} textAnchor="middle" fill="#78928b" fontSize={11} fontFamily="DM Mono, monospace">
          {tick === 0 ? '0' : signed(tick, 2)}
        </text>
      </g>)}
      <text x={x(0)} y={22} textAnchor="middle" fill="#8da29d" fontSize={11}>no effect</text>

      {/* Primary: band + point estimate */}
      <rect x={x(low)} y={50} width={x(high) - x(low)} height={16} rx={4} fill="#6bf5bb" opacity={0.22}/>
      <line x1={x(low)} x2={x(low)} y1={46} y2={70} stroke="#6bf5bb" strokeWidth={2}/>
      <line x1={x(high)} x2={x(high)} y1={46} y2={70} stroke="#6bf5bb" strokeWidth={2}/>
      <circle cx={x(primary.median_shift_hodges_lehmann)} cy={58} r={6} fill="#6bf5bb" stroke="#0b1815" strokeWidth={2}/>
      <text x={52} y={42} fill="#edf8f4" fontSize={12} fontWeight={600}>Traditional antidiabetic use</text>
      <text x={x(high) + 14} y={63} fill="#8da29d" fontSize={11} fontFamily="DM Mono, monospace">
        {signed(primary.median_shift_hodges_lehmann, 4)} · p {primary.p_value.toFixed(2)}
      </text>

      {/* Negative control: point estimate only */}
      <circle cx={x(control.median_shift_hodges_lehmann)} cy={110} r={6} fill="#f0ba69" stroke="#0b1815" strokeWidth={2}/>
      <text x={52} y={96} fill="#edf8f4" fontSize={12} fontWeight={600}>Negative control (toothache)</text>
      <text x={x(control.median_shift_hodges_lehmann) + 16} y={115} fill="#8da29d" fontSize={11} fontFamily="DM Mono, monospace">
        {signed(control.median_shift_hodges_lehmann, 4)} · p {control.p_value.toFixed(2)}
      </text>

      <text x={350} y={180} textAnchor="middle" fill="#78928b" fontSize={11}>
        Median shift in predicted pChEMBL (plant-level Hodges–Lehmann)
      </text>
    </svg>
    <figcaption>
      The shaded band is the prespecified family-clustered 95% interval; it spans zero, and it is
      narrow — so this bounds the effect rather than merely failing to find one. The negative
      control, a traditional use with no expected PTP1B relevance, sits <strong>further</strong> from
      zero than the hypothesis under test. No bootstrap interval was prespecified for the control,
      so it is shown as a point estimate only.
    </figcaption>
  </figure>
}

function PlantExplorer() {
  const [query, setQuery] = useState('')
  const [selected, setSelected] = useState<EthnoPlant>(
    () => summary.plants.find((plant) => plant.species === DEFAULT_SPECIES) ?? summary.plants[0],
  )

  const matches = useMemo(() => {
    const needle = query.trim().toLowerCase()
    if (!needle) return []
    return summary.plants
      .filter((plant) => plant.species.includes(needle) || plant.family.includes(needle))
      .slice(0, 8)
  }, [query])

  const cut = summary.high_potential_definition.predicted_pchembl_at_least
  const outOfDomain = selected.compounds - selected.inDomain

  return <section className="panel plant-explorer">
    <div className="panel-head">
      <div>
        <span className="section-kicker">PER-PLANT READOUT</span>
        <h2>What the model can and cannot say about one plant</h2>
      </div>
      <Leaf size={19}/>
    </div>

    <div className="plant-search">
      <label htmlFor="plant-query"><Search size={15}/> Search {summary.totals.plants.toLocaleString()} analysed plants</label>
      <div>
        <input
          id="plant-query"
          type="search"
          value={query}
          placeholder="e.g. coix, eucalyptus, rosaceae"
          onChange={(event) => setQuery(event.target.value)}
          autoComplete="off"
        />
        {query && <button type="button" className="plant-clear" onClick={() => setQuery('')} aria-label="Clear search"><X size={14}/></button>}
      </div>
      {query.trim() && <ul className="plant-matches">
        {matches.length === 0 && <li className="plant-empty">No analysed plant matches that name.</li>}
        {matches.map((plant) => <li key={plant.species}>
          <button type="button" onClick={() => { setSelected(plant); setQuery('') }}>
            <strong>{titleCase(plant.species)}</strong>
            <span>{titleCase(plant.family)}</span>
            {plant.isCase && <i className="ethno-badge case">traditional antidiabetic</i>}
          </button>
        </li>)}
      </ul>}
    </div>

    <div className="plant-result">
      <header>
        <div>
          <h3>{titleCase(selected.species)}</h3>
          <span>{titleCase(selected.family)}</span>
        </div>
        {selected.isCase
          ? <i className="ethno-badge case">traditional antidiabetic use</i>
          : <i className="ethno-badge control">no recorded antidiabetic use</i>}
      </header>

      <div className="domain-bar" role="img" aria-label={
        `${selected.compounds} compounds scored: ${selected.inDomain} inside the applicability domain, `
        + `${outOfDomain} outside, ${selected.highPotential} high potential.`
      }>
        {selected.inDomain > 0 && <span
          className="seg in"
          style={{ flexGrow: selected.inDomain }}
        >{selected.inDomain} in-domain</span>}
        {outOfDomain > 0 && <span
          className="seg out"
          style={{ flexGrow: outOfDomain }}
        >{outOfDomain} out-of-domain</span>}
      </div>

      <p className="plant-verdict">
        Of the <strong>{selected.compounds}</strong> natural products NPASS reports for this plant,
        {' '}<strong>{selected.inDomain}</strong> fall inside the model&apos;s applicability domain
        {selected.highPotential > 0
          ? <> and <strong>{selected.highPotential}</strong> of those reach predicted pChEMBL ≥ {cut}</>
          : <> and <strong>none</strong> reach predicted pChEMBL ≥ {cut}</>}.
        {' '}The remaining <strong>{outOfDomain}</strong> are chemistry this model was not trained on,
        so their scores are not evidence either way.
        {' '}<em>Nothing here is clinical proof, and nothing here is measured activity.</em>
      </p>

      {selected.activities && <dl className="plant-provenance">
        <div><dt>Recorded traditional use</dt><dd>{selected.activities}</dd></div>
        {selected.countries && <div><dt>Country</dt><dd>{selected.countries}</dd></div>}
        {selected.references && <div><dt>Duke reference</dt><dd>{selected.references}</dd></div>}
      </dl>}

      <small className="plant-footnote">
        &ldquo;High potential&rdquo; = predicted pChEMBL ≥ {cut} (1 µM) <strong>and</strong> inside the
        applicability domain. Across all {summary.totals.plants.toLocaleString()} plants, only
        {' '}{summary.totals.compounds_high_potential} of {summary.totals.compounds_scored.toLocaleString()} scored
        compounds clear both bars. {summary.high_potential_definition.note}
      </small>
    </div>
  </section>
}

export function ScientificValidationWorkspace() {
  const primary = results.primary_outcome
  const conditions = [
    {
      label: `Mann-Whitney p < 0.05`,
      detail: `observed p = ${primary.p_value.toFixed(4)}`,
      met: results.support_rule.conditions.p_below_alpha,
    },
    {
      label: 'Family-clustered 95% CI excludes 0',
      detail: `observed ${signed(primary.family_clustered_ci_95[0], 4)} to ${signed(primary.family_clustered_ci_95[1], 4)}`,
      met: results.support_rule.conditions.ci_excludes_zero,
    },
    {
      label: 'Negative-control shift smaller than primary',
      detail: `control ${signed(results.negative_control.median_shift_hodges_lehmann, 4)} vs primary ${signed(primary.median_shift_hodges_lehmann, 4)}`,
      met: results.support_rule.conditions.negative_control_smaller,
    },
  ]

  const inDomainShare = summary.totals.compounds_in_domain / summary.totals.compounds_scored

  return <>
    <header className="page-title experiment-title">
      <div>
        <span className="section-kicker">PRESPECIFIED ETHNOPHARMACOLOGICAL TEST · {results.protocol_version}</span>
        <h1>Does traditional use predict measurable activity?</h1>
        <p>{results.question}</p>
      </div>
      <div className="experiment-verdict">
        <AlertTriangle/>
        <div><strong>{results.support_rule.verdict}</strong><span>0 of 3 prespecified conditions met</span></div>
      </div>
    </header>

    <section className="experiment-context-strip">
      <div><span>PLANTS ANALYSED</span><strong>{results.cohort.analysed_plants.toLocaleString()}</strong><small>{results.cohort.families} botanical families</small></div>
      <div><span>COMPOUNDS SCORED</span><strong>{results.cohort.compounds_scored.toLocaleString()}</strong><small>{results.cohort.compounds_unusable} rejected, never imputed</small></div>
      <div><span>TRADITIONAL ANTIDIABETIC</span><strong>{primary.case_plants}</strong><small>vs {primary.control_plants.toLocaleString()} Duke-covered controls</small></div>
      <div><span>INSIDE APPLICABILITY DOMAIN</span><strong>{pct(summary.totals.compounds_in_domain, summary.totals.compounds_scored)}</strong><small>{summary.totals.compounds_in_domain.toLocaleString()} of {summary.totals.compounds_scored.toLocaleString()} compounds</small></div>
    </section>

    <section className="panel convergence-timeline">
      <div className="panel-head">
        <div>
          <span className="section-kicker">WHY THIS IS A TEST AND NOT A RETROFIT</span>
          <h2>The protocol was public before the outcome existed</h2>
        </div>
        <span className="artifact-badge">VERIFIABLE</span>
      </div>
      <ol>
        {timeline.map((entry) => <li key={entry.sha} className={entry.isOutcome ? 'outcome' : undefined}>
          <GitCommitHorizontal size={16}/>
          <time>{entry.time}</time>
          <a href={`${REPO}/${entry.sha}`} target="_blank" rel="noreferrer">
            <code>{entry.sha.slice(0, 7)}</code><ExternalLink size={12}/>
          </a>
          <p>{entry.what}</p>
        </li>)}
      </ol>
      <p className="timeline-note">
        Two hours and ten minutes separate the frozen protocol from the first number. Open the
        first commit and read the support rule, the α, the seed and the negative-control terms as
        they stood before any score was computed — then compare them with the verdict below.
      </p>
    </section>

    <div className="experiment-grid">
      <section className="panel experiment-chart-panel">
        <div className="panel-head">
          <div><span className="section-kicker">PRIMARY OUTCOME</span><h2>Effect size with its prespecified interval</h2></div>
        </div>
        <IntervalPlot/>
      </section>

      <aside className="panel metric-inspector">
        <span className="section-kicker">PRESPECIFIED SUPPORT RULE</span>
        <p className="rule-intro">All three conditions were required. None were met.</p>
        <ul className="condition-list">
          {conditions.map((condition) => <li key={condition.label} className={condition.met ? 'met' : 'unmet'}>
            <span aria-hidden="true">{condition.met ? '✓' : '✕'}</span>
            <div>
              <strong>{condition.label}</strong>
              <small>{condition.detail}</small>
            </div>
            <i className="sr-only">{condition.met ? 'met' : 'not met'}</i>
          </li>)}
        </ul>
        <div className="difference-callout inconclusive">
          <span>VERDICT, AS WRITTEN IN ADVANCE</span>
          <strong>{results.support_rule.verdict}</strong>
          <small>Reported unchanged. No secondary analysis was promoted, and the negative control was not dropped.</small>
        </div>
      </aside>
    </div>

    <section className="panel caveat-panel">
      <div className="panel-head">
        <div><span className="section-kicker">THE BOUND ON OUR OWN CLAIM</span><h2>The caveat we do not bury</h2></div>
        <ShieldAlert size={19}/>
      </div>
      <div className="caveat-layout">
        <div className="negative-result">
          <AlertTriangle/>
          <div>
            <strong>{pct(summary.totals.compounds_scored - summary.totals.compounds_in_domain, summary.totals.compounds_scored)} of this chemistry is outside the model&apos;s reach</strong>
            <p>
              A model scoring chemistry it was never trained on regresses toward its own mean, and a
              flattened score cannot express a difference even where one exists. So the honest claim
              is bounded: traditional antidiabetic use does not predict PTP1B activity
              {' '}<em>as this model estimates it</em>. That is not evidence that traditional
              knowledge lacks pharmacological content.
            </p>
          </div>
        </div>
        <div className="sensitivity-card">
          <span className="section-kicker">PRESPECIFIED SENSITIVITY · {sensitivity.prespecified_as}</span>
          <h3>Restricting to in-domain chemistry does not rescue it</h3>
          <dl>
            <div><dt>Plants retained</dt><dd>{sensitivity.analysed_plants.toLocaleString()}</dd></div>
            <div><dt>Median shift</dt><dd>{signed(sensitivity.primary_outcome.median_shift_hodges_lehmann, 4)}</dd></div>
            <div><dt>Mann-Whitney p</dt><dd>{sensitivity.primary_outcome.p_value.toFixed(4)}</dd></div>
            <div><dt>95% CI</dt><dd>{signed(sensitivity.primary_outcome.family_clustered_ci_95[0], 4)} to {signed(sensitivity.primary_outcome.family_clustered_ci_95[1], 4)}</dd></div>
            <div><dt>Negative control</dt><dd>{signed(sensitivity.negative_control.median_shift_hodges_lehmann, 4)}</dd></div>
          </dl>
          <p>
            The shift turns positive but stays well inside noise, and the negative control rises with
            it — still the larger of the two. Narrowing to the {pct(summary.totals.compounds_in_domain, summary.totals.compounds_scored)} the
            model actually covers reproduces the null on a smaller cohort.
            {!sensitivity.convergence_supported && ' Convergence remains unsupported.'}
          </p>
        </div>
      </div>
    </section>

    <PlantExplorer/>

    <section className="panel limitations-panel">
      <div className="panel-head">
        <div><span className="section-kicker">DECLARED IN THE PROTOCOL, BEFORE THE RESULT</span><h2>Limitations</h2></div>
        <span className="artifact-badge">{Math.round(inDomainShare * 100)}% IN DOMAIN</span>
      </div>
      <ul className="limitations-list">
        {results.declared_limitations.map((limitation) => <li key={limitation}>{limitation}</li>)}
      </ul>
      <footer className="integrity-footer">
        <div>
          <strong>Traditional knowledge handling</strong>
          <p>
            Every traditional-use record keeps its original Duke reference and country, shown above
            for any plant that has one. This analysis reports convergence with a model prediction;
            it never states that traditional knowledge is validated or invalidated by one. No
            genetic resources were accessed — only published literature records.
          </p>
        </div>
        <div>
          <strong>Input checksums</strong>
          <ul className="checksum-list">
            {Object.entries(results.inputs_sha256).map(([name, digest]) => <li key={name}>
              <span>{name}</span><code>{digest.slice(0, 16)}…</code>
            </li>)}
          </ul>
        </div>
      </footer>
    </section>
  </>
}
