import { afterEach, describe, expect, it, vi } from 'vitest'
import { predictMolecule } from './api'

describe('prediction API service', () => {
  afterEach(() => vi.unstubAllGlobals())

  it('sends a typed JSON request to the real endpoint contract', async () => {
    const fetchMock = vi.fn().mockResolvedValue(new Response(JSON.stringify({ canonical_smiles: 'CCO' }), {
      status: 200,
      headers: { 'Content-Type': 'application/json' },
    }))
    vi.stubGlobal('fetch', fetchMock)

    await predictMolecule('CCO')

    expect(fetchMock).toHaveBeenCalledWith('/api/predict', expect.objectContaining({
      method: 'POST',
      body: JSON.stringify({ smiles: 'CCO' }),
    }))
  })

  it('preserves invalid-SMILES status and error code', async () => {
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue(new Response(JSON.stringify({
      detail: { code: 'invalid_smiles', message: 'The supplied SMILES is not a valid molecule.' },
    }), { status: 422, headers: { 'Content-Type': 'application/json' } })))

    await expect(predictMolecule('C1CC')).rejects.toEqual(expect.objectContaining({
      status: 422,
      code: 'invalid_smiles',
    }))
  })

  it('reports network failures without inventing a result', async () => {
    vi.stubGlobal('fetch', vi.fn().mockRejectedValue(new TypeError('offline')))
    await expect(predictMolecule('CCO')).rejects.toMatchObject({ code: 'network_error', status: 0 })
  })

  it('turns a proxy 502 into a useful unavailable-backend error', async () => {
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue(new Response('Bad Gateway', { status: 502 })))
    await expect(predictMolecule('CCO')).rejects.toMatchObject({
      code: 'backend_unavailable',
      status: 502,
      message: 'The PhytoTruth API is unavailable. Confirm that the FastAPI service is running.',
    })
  })
})
