import { copyFileSync, existsSync, mkdirSync } from 'node:fs'
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const frontend = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const artifacts = [
  ['evidence_reversal_chembl335', 'evidence-reversal-results.json'],
  ['robust_multitarget_evaluation', 'robust-multitarget-results.json'],
  ['evidence_intelligence', 'evidence-intelligence-results.json'],
]

for (const [experiment, filename] of artifacts) {
  const source = resolve(frontend, '..', 'experiments', experiment, 'results.json')
  const destination = resolve(frontend, 'src', 'data', filename)
  if (!existsSync(source)) {
    throw new Error(`Missing executed scientific artifact: ${source}`)
  }
  mkdirSync(dirname(destination), { recursive: true })
  copyFileSync(source, destination)
  console.log(`Synchronized executed experiment artifact: ${destination}`)
}
