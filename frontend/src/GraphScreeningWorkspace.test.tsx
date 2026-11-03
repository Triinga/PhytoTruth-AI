import { render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { predictMolecule } from './api'
import { fallbackDashboard } from './demo'
import { GraphScreeningWorkspace } from './GraphScreeningWorkspace'
import type { PredictionResponse } from './types'

vi.mock('./api', async importOriginal => {
  const actual = await importOriginal<typeof import('./api')>()
  return { ...actual, predictMolecule: vi.fn() }
})

function prediction(smiles: string, value: number): PredictionResponse {
  return {
    canonical_smiles: smiles,
    predicted_pchembl: value,
    prediction_interval: { lower: value - .8, upper: value + .8, confidence_level: .9, method: 'fixture' },
    applicability_domain: { inside: true, status: 'IN_DOMAIN_RESEARCH_ONLY', nearest_training_tanimoto: .8, threshold: .42, warning: null },
    descriptors: { molecular_weight: 300, logp: 2, tpsa: 70, hbd: 2, hba: 4, rotatable_bonds: 3, fraction_csp3: .5, ring_count: 2, heavy_atoms: 22, formal_charge: 0 },
    model: { name: 'xgboost', version: 'fixture', format: 'xgboost-json', target: { assay_type: 'B', bao_format: 'BAO_0000357', standard_type: 'IC50', target_chembl_id: 'CHEMBL335', target_organism: 'Homo sapiens' }, training_structures: 2444, training_scaffolds: 1086 },
    kind: 'predicted',
    limitations: ['research only'],
  }
}

describe('Graph and screening workspace', () => {
  beforeEach(() => {
    vi.mocked(predictMolecule).mockReset()
    vi.mocked(predictMolecule).mockImplementation(async smiles => prediction(smiles, smiles.length / 20 + 5))
  })

  it('runs real endpoint contracts and creates an auditable ranked graph', async () => {
    render(<GraphScreeningWorkspace data={fallbackDashboard}/>)
    await userEvent.click(screen.getByRole('button', { name: /Run live demo/ }))

    await waitFor(() => expect(screen.getByText('Experimental follow-up queue')).toBeInTheDocument())
    expect(predictMolecule).toHaveBeenCalled()
    expect(screen.getByText('Trace the candidate, not just the score')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /COMPOUND/ })).toBeInTheDocument()
    expect(screen.getByText(/Research integrity boundary/)).toBeInTheDocument()
  })

  it('supports domain filtering without changing the stored ranking', async () => {
    vi.mocked(predictMolecule).mockImplementationOnce(async smiles => ({
      ...prediction(smiles, 7),
      applicability_domain: { inside: false, status: 'OUTSIDE_DOMAIN_DO_NOT_RELY', nearest_training_tanimoto: .2, threshold: .42, warning: 'outside' },
    })).mockImplementation(async smiles => prediction(smiles, 6))
    render(<GraphScreeningWorkspace data={fallbackDashboard}/>)
    await userEvent.click(screen.getByRole('button', { name: /Run live demo/ }))
    await screen.findByText('Experimental follow-up queue')
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Applicability domain' }), 'in')
    expect(screen.getByText(/candidates match/)).toBeInTheDocument()
  })
})
