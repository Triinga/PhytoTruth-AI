import { fireEvent, render, screen, waitFor, within } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { fetchMoleculeSvg, predictMolecule } from './api'
import { fallbackDashboard } from './demo'
import { PredictionStructure, PredictionWorkspace } from './PredictionWorkspace'
import type { PredictionResponse } from './types'

vi.mock('./api', async importOriginal => {
  const actual = await importOriginal<typeof import('./api')>()
  return { ...actual, predictMolecule: vi.fn(), fetchMoleculeSvg: vi.fn() }
})

const rdkitSvg = `<?xml version="1.0" encoding="iso-8859-1"?>
<svg xmlns="http://www.w3.org/2000/svg" width="520px" height="260px" viewBox="0 0 520 260">
  <path d="M 20,130 L 500,130" style="fill:none;stroke:#000000;stroke-width:2px"/>
  <text x="245" y="120" style="fill:#000000">O</text>
</svg>`

// This is the recorded CCO response used to test rendering; production values
// are always requested from FastAPI and never sourced from this fixture.
const response: PredictionResponse = {
  canonical_smiles: 'CCO',
  predicted_pchembl: 4.849597454071045,
  prediction_interval: { lower: 2.988195328965653, upper: 6.710999579176437, confidence_level: .9, method: 'cross-validated absolute-residual conformal heuristic' },
  applicability_domain: { inside: false, status: 'OUTSIDE_DOMAIN_DO_NOT_RELY', nearest_training_tanimoto: .21052631578947367, threshold: .4225352112676056, warning: 'Outside applicability domain; interval calibration is not guaranteed.' },
  descriptors: { molecular_weight: 46.069, logp: -.001, tpsa: 20.23, hbd: 1, hba: 1, rotatable_bonds: 0, fraction_csp3: 1, ring_count: 0, heavy_atoms: 3, formal_charge: 0 },
  model: { name: 'xgboost', version: 'chembl335-xgboost-c20d771793ac', format: 'xgboost-json', target: { assay_type: 'B', bao_format: 'BAO_0000357', standard_type: 'IC50', target_chembl_id: 'CHEMBL335', target_organism: 'Homo sapiens' }, training_structures: 2444, training_scaffolds: 1086 },
  kind: 'predicted',
  limitations: ['Research prioritization only.'],
}

describe('PredictionWorkspace', () => {
  beforeEach(() => {
    vi.mocked(predictMolecule).mockReset()
    vi.mocked(fetchMoleculeSvg).mockReset().mockResolvedValue(rdkitSvg)
  })

  it('validates blank input before calling the API', async () => {
    render(<PredictionWorkspace candidates={fallbackDashboard.candidates}/>)
    await userEvent.click(screen.getByRole('button', { name: /run prediction/i }))
    expect(screen.getByRole('alert')).toHaveTextContent('Enter a SMILES')
    expect(predictMolecule).not.toHaveBeenCalled()
  })

  it('renders an out-of-domain model response transparently', async () => {
    vi.mocked(predictMolecule).mockResolvedValue(response)
    render(<PredictionWorkspace candidates={fallbackDashboard.candidates}/>)
    await userEvent.type(screen.getByLabelText('Molecular SMILES'), 'CCO')
    await userEvent.click(screen.getByRole('button', { name: /run prediction/i }))

    expect(await screen.findByText('4.850')).toBeInTheDocument()
    expect(screen.getByText('Outside modeled chemistry')).toBeInTheDocument()
    expect(screen.getByText(/Do not treat this prediction as reliable/)).toBeInTheDocument()
    expect(screen.getByText('chembl335-xgboost-c20d771793ac')).toBeInTheDocument()
    expect(predictMolecule).toHaveBeenCalledWith('CCO', expect.any(AbortSignal))
  })

  it('invalidates a displayed prediction as soon as the SMILES changes', async () => {
    vi.mocked(predictMolecule).mockResolvedValue(response)
    render(<PredictionWorkspace candidates={fallbackDashboard.candidates}/>)
    const input = screen.getByLabelText('Molecular SMILES')
    await userEvent.type(input, 'CCO')
    await userEvent.click(screen.getByRole('button', { name: /run prediction/i }))
    expect(await screen.findByText('4.850')).toBeInTheDocument()

    await userEvent.type(input, 'N')

    expect(screen.queryByText('4.850')).not.toBeInTheDocument()
    expect(screen.getByText('No prediction yet')).toBeInTheDocument()
  })

  it('prevents an older asynchronous response from replacing a newer result', async () => {
    let resolveFirst!: (value: PredictionResponse) => void
    let resolveSecond!: (value: PredictionResponse) => void
    const first = new Promise<PredictionResponse>(resolve => { resolveFirst = resolve })
    const second = new Promise<PredictionResponse>(resolve => { resolveSecond = resolve })
    vi.mocked(predictMolecule).mockReturnValueOnce(first).mockReturnValueOnce(second)
    render(<PredictionWorkspace candidates={fallbackDashboard.candidates}/>)
    const input = screen.getByLabelText('Molecular SMILES')

    await userEvent.type(input, 'CCO')
    await userEvent.click(screen.getByRole('button', { name: /run prediction/i }))
    await userEvent.clear(input)
    await userEvent.type(input, 'CCN')
    await userEvent.click(screen.getByRole('button', { name: /run prediction/i }))

    resolveSecond({ ...response, canonical_smiles: 'CCN', predicted_pchembl: 6.123 })
    expect(await screen.findByText('6.123')).toBeInTheDocument()
    resolveFirst(response)
    await Promise.resolve()

    expect(screen.getByText('6.123')).toBeInTheDocument()
    expect(screen.queryByText('4.850')).not.toBeInTheDocument()
  })
})

describe('PredictionStructure', () => {
  beforeEach(() => vi.mocked(fetchMoleculeSvg).mockReset().mockResolvedValue(rdkitSvg))

  it('renders RDKit markup as inline vector SVG and never as raster content', async () => {
    const { container } = render(<PredictionStructure result={response}/>)
    await waitFor(() => expect(container.querySelector('svg.rdkit-vector-svg')).not.toBeNull())
    const svg = container.querySelector('svg.rdkit-vector-svg')

    expect(svg).toHaveAttribute('viewBox', '97.5 48.75 325 162.5')
    expect(svg).toHaveAttribute('preserveAspectRatio', 'xMidYMid meet')
    expect(svg).toHaveAttribute('width', '100%')
    expect(container.querySelector('img,canvas')).not.toBeInTheDocument()
    expect(svg?.querySelector('path')?.getAttribute('style')).toContain('stroke:#dceee8')
  })

  it('zooms the vector viewBox and restores the fitted view', async () => {
    render(<PredictionStructure result={response}/>)
    const zoomIn = await screen.findByRole('button', { name: 'Zoom in molecule' })
    const viewport = screen.getByTestId('molecule-viewport')
    const svg = viewport.querySelector('svg')!

    expect(svg).toHaveAttribute('viewBox', '97.5 48.75 325 162.5')
    await userEvent.click(zoomIn)
    expect(svg).toHaveAttribute('viewBox', '130 65 260 130')
    expect(screen.getByText('200%')).toBeInTheDocument()

    await userEvent.click(screen.getByRole('button', { name: 'Fit molecule to view' }))
    expect(svg).toHaveAttribute('viewBox', '0 0 520 260')
    expect(screen.getByText('100%')).toBeInTheDocument()
  })

  it('supports cursor-centered wheel zoom and drag panning inside the viewer', async () => {
    render(<PredictionStructure result={response}/>)
    const viewport = await screen.findByTestId('molecule-viewport')
    const svg = viewport.querySelector('svg')!
    vi.spyOn(viewport, 'getBoundingClientRect').mockReturnValue({
      x: 0, y: 0, left: 0, top: 0, right: 400, bottom: 200, width: 400, height: 200, toJSON: () => ({}),
    })

    fireEvent.wheel(viewport, { deltaY: -100, clientX: 100, clientY: 50 })
    expect(svg).toHaveAttribute('viewBox', '113.75 56.875 260 130')

    fireEvent.pointerDown(viewport, { button: 0, pointerId: 1, clientX: 10, clientY: 15 })
    fireEvent.pointerMove(viewport, { pointerId: 1, clientX: 50, clientY: 35 })
    fireEvent.pointerUp(viewport, { pointerId: 1 })
    expect(svg).toHaveAttribute('viewBox', '87.75 43.875 260 130')
  })

  it('automatically fits a newly loaded molecule', async () => {
    const { rerender } = render(<PredictionStructure result={response}/>)
    await userEvent.click(await screen.findByRole('button', { name: 'Zoom in molecule' }))
    expect(screen.getByText('200%')).toBeInTheDocument()

    rerender(<PredictionStructure result={{ ...response, canonical_smiles: 'CCN' }}/>)

    expect(await screen.findByText('160%')).toBeInTheDocument()
    expect(vi.mocked(fetchMoleculeSvg)).toHaveBeenLastCalledWith('CCN', expect.any(AbortSignal))
  })

  it('opens a near-fullscreen viewer with reset, close, and Escape handling', async () => {
    render(<PredictionStructure result={response}/>)
    await userEvent.click(await screen.findByRole('button', { name: 'Open expanded molecule viewer' }))
    const dialog = screen.getByRole('dialog', { name: /expanded molecular structure/i })
    expect(document.body.style.overflow).toBe('hidden')
    expect(within(dialog).getByRole('button', { name: 'Reset molecule view' })).toBeInTheDocument()

    await userEvent.click(within(dialog).getByRole('button', { name: 'Zoom in molecule' }))
    expect(within(dialog).getByText('125%')).toBeInTheDocument()
    await userEvent.click(within(dialog).getByRole('button', { name: 'Reset molecule view' }))
    expect(within(dialog).getByText('100%')).toBeInTheDocument()

    fireEvent.keyDown(window, { key: 'Escape' })
    await waitFor(() => expect(screen.queryByRole('dialog')).not.toBeInTheDocument())
    expect(document.body.style.overflow).toBe('')

    await userEvent.click(screen.getByRole('button', { name: 'Open expanded molecule viewer' }))
    await userEvent.click(screen.getByRole('button', { name: 'Close expanded molecule viewer' }))
    expect(screen.queryByRole('dialog')).not.toBeInTheDocument()
  })
})
