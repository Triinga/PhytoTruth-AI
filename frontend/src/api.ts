import type { Dashboard, PredictionResponse, ScreeningArtifactResponse } from './types'

export const API_BASE = import.meta.env.VITE_API_URL ?? ''

export class ApiError extends Error {
  constructor(
    message: string,
    readonly status: number,
    readonly code?: string,
  ) {
    super(message)
    this.name = 'ApiError'
  }
}

export function apiUrl(path: string): string {
  return `${API_BASE}${path}`
}

export function apiDocsUrl(): string {
  return API_BASE ? `${API_BASE}/docs` : `${window.location.protocol}//${window.location.hostname}:8000/docs`
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  let response: Response
  try {
    response = await fetch(apiUrl(path), init)
  } catch {
    throw new ApiError('The PhytoTruth API could not be reached. Confirm that the FastAPI service is running.', 0, 'network_error')
  }

  if (!response.ok) {
    let message = `The API returned HTTP ${response.status}.`
    let code: string | undefined
    let hasServerDetail = false
    try {
      const payload = await response.json() as { detail?: string | { message?: string; code?: string } }
      if (typeof payload.detail === 'string') {
        message = payload.detail
        hasServerDetail = true
      }
      else if (payload.detail) {
        message = payload.detail.message ?? message
        code = payload.detail.code
        hasServerDetail = Boolean(payload.detail.message || payload.detail.code)
      }
    } catch {
      // Preserve the status-based error if the server did not return JSON.
    }
    if (response.status >= 500 && !hasServerDetail) {
      message = 'The PhytoTruth API is unavailable. Confirm that the FastAPI service is running.'
      code = 'backend_unavailable'
    }
    throw new ApiError(message, response.status, code)
  }
  return response.json() as Promise<T>
}

export function getDashboard(signal?: AbortSignal): Promise<Dashboard> {
  return request<Dashboard>('/api/dashboard', { signal })
}

export function predictMolecule(smiles: string, signal?: AbortSignal): Promise<PredictionResponse> {
  return request<PredictionResponse>('/api/predict', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ smiles }),
    signal,
  })
}

export function getScreeningCandidates(limit = 100, signal?: AbortSignal): Promise<ScreeningArtifactResponse> {
  return request<ScreeningArtifactResponse>(`/api/screening/candidates?limit=${encodeURIComponent(limit)}`, { signal })
}

export function moleculeSvgUrl(smiles: string): string {
  return apiUrl(`/api/molecule/svg?smiles=${encodeURIComponent(smiles)}`)
}

export async function fetchMoleculeSvg(smiles: string, signal?: AbortSignal): Promise<string> {
  let response: Response
  try {
    response = await fetch(moleculeSvgUrl(smiles), { signal })
  } catch (cause) {
    if (cause instanceof DOMException && cause.name === 'AbortError') throw cause
    throw new ApiError('The RDKit structure renderer could not be reached.', 0, 'renderer_unavailable')
  }
  if (!response.ok) {
    throw new ApiError(
      response.status === 422 ? 'RDKit could not render this molecular structure.' : 'The RDKit structure renderer is unavailable.',
      response.status,
      response.status === 422 ? 'invalid_smiles' : 'renderer_unavailable',
    )
  }
  const markup = await response.text()
  if (!markup.trim().startsWith('<?xml') && !markup.trim().startsWith('<svg')) {
    throw new ApiError('The structure renderer returned invalid SVG.', 502, 'invalid_svg')
  }
  return markup
}
