import { copyFileSync, createReadStream, existsSync, mkdirSync, writeFileSync } from 'node:fs'
import { createInterface } from 'node:readline'
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const frontend = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const repo = resolve(frontend, '..')
const dataDir = resolve(frontend, 'src', 'data')

const artifacts = [
  ['evidence_reversal_chembl335', 'results.json', 'evidence-reversal-results.json'],
  ['robust_multitarget_evaluation', 'results.json', 'robust-multitarget-results.json'],
  ['evidence_intelligence', 'results.json', 'evidence-intelligence-results.json'],
  ['ethno_convergence', 'results.json', 'ethno-convergence-results.json'],
  ['ethno_convergence', 'sensitivity_in_domain.json', 'ethno-convergence-sensitivity.json'],
]

for (const [experiment, file, filename] of artifacts) {
  const source = resolve(repo, 'experiments', experiment, file)
  const destination = resolve(dataDir, filename)
  if (!existsSync(source)) {
    throw new Error(`Missing executed scientific artifact: ${source}`)
  }
  mkdirSync(dirname(destination), { recursive: true })
  copyFileSync(source, destination)
  console.log(`Synchronized executed experiment artifact: ${destination}`)
}

// ---------------------------------------------------------------------------
// Per-plant summary derived from the executed convergence artifacts.
//
// "High potential" is a prespecified, documented cut applied to model output:
// predicted pChEMBL >= 6.0 (1 uM) AND inside the model applicability domain.
// Both bars must be cleared, because a high score on chemistry the model was
// not trained on is not evidence. This is a research-prioritisation label, not
// a claim of measured activity.
// ---------------------------------------------------------------------------
const HIGH_POTENTIAL_PCHEMBL = 6
const plantCsv = resolve(repo, 'experiments', 'ethno_convergence', 'plant_level_scores.csv')
const compoundCsv = resolve(repo, 'experiments', 'ethno_convergence', 'compound_level_scores.csv')

for (const path of [plantCsv, compoundCsv]) {
  if (!existsSync(path)) {
    throw new Error(`Missing executed convergence table: ${path}`)
  }
}

/** Minimal RFC 4180 reader: handles quoted fields containing commas and "" escapes. */
function parseCsvLine(line) {
  const out = []
  let field = ''
  let quoted = false
  for (let i = 0; i < line.length; i += 1) {
    const char = line[i]
    if (quoted) {
      if (char === '"') {
        if (line[i + 1] === '"') { field += '"'; i += 1 } else { quoted = false }
      } else { field += char }
    } else if (char === '"') { quoted = true } else if (char === ',') { out.push(field); field = '' } else { field += char }
  }
  out.push(field)
  return out
}

async function readCsv(path, onRow) {
  const reader = createInterface({ input: createReadStream(path, 'utf8'), crlfDelay: Infinity })
  let header = null
  for await (const line of reader) {
    if (!line) continue
    const cells = parseCsvLine(line)
    if (!header) { header = cells; continue }
    const row = {}
    header.forEach((key, index) => { row[key] = cells[index] ?? '' })
    onRow(row)
  }
}

const highPotentialByPlant = new Map()
let compoundsScored = 0
let compoundsInDomain = 0
let compoundsHighPotential = 0

await readCsv(compoundCsv, (row) => {
  compoundsScored += 1
  const inDomain = row.inside_applicability_domain === 'True'
  if (inDomain) compoundsInDomain += 1
  const value = Number.parseFloat(row.predicted_pchembl)
  if (inDomain && Number.isFinite(value) && value >= HIGH_POTENTIAL_PCHEMBL) {
    compoundsHighPotential += 1
    highPotentialByPlant.set(row.species, (highPotentialByPlant.get(row.species) ?? 0) + 1)
  }
})

const number = (value) => {
  const parsed = Number.parseFloat(value)
  return Number.isFinite(parsed) ? Math.round(parsed * 1e4) / 1e4 : null
}

const plants = []
await readCsv(plantCsv, (row) => {
  const plant = {
    species: row.species,
    family: row.family,
    isCase: row.is_case === 'True',
    isNegativeControl: row.is_negative_control === 'True',
    compounds: Number.parseInt(row.n_compounds_scored, 10) || 0,
    inDomain: Number.parseInt(row.n_in_domain, 10) || 0,
    highPotential: highPotentialByPlant.get(row.species) ?? 0,
    score: number(row.score),
    scoreInDomain: number(row.score_in_domain),
  }
  // Attribution is only carried for plants that actually have a recorded use.
  if (row.traditional_activities) plant.activities = row.traditional_activities
  if (row.duke_references) plant.references = row.duke_references
  if (row.duke_countries) plant.countries = row.duke_countries
  plants.push(plant)
})

plants.sort((left, right) => left.species.localeCompare(right.species))

const summary = {
  generated_from: {
    plant_level_scores: 'experiments/ethno_convergence/plant_level_scores.csv',
    compound_level_scores: 'experiments/ethno_convergence/compound_level_scores.csv',
  },
  high_potential_definition: {
    predicted_pchembl_at_least: HIGH_POTENTIAL_PCHEMBL,
    requires_inside_applicability_domain: true,
    note: 'Research-prioritisation label applied to model output. Not measured activity.',
  },
  totals: {
    plants: plants.length,
    families: new Set(plants.map((plant) => plant.family)).size,
    case_plants: plants.filter((plant) => plant.isCase).length,
    compounds_scored: compoundsScored,
    compounds_in_domain: compoundsInDomain,
    compounds_high_potential: compoundsHighPotential,
  },
  plants,
}

const summaryPath = resolve(dataDir, 'ethno-plants.json')
writeFileSync(summaryPath, `${JSON.stringify(summary)}\n`, 'utf8')
console.log(
  `Derived per-plant convergence summary: ${summaryPath} `
  + `(${summary.totals.plants} plants, ${compoundsHighPotential} high-potential compounds)`,
)
