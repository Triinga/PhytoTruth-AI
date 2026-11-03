import { fireEvent, render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { getDashboard } from './api'
import App from './App'
import { fallbackDashboard } from './demo'

vi.mock('./api', async importOriginal => {
  const actual = await importOriginal<typeof import('./api')>()
  return { ...actual, getDashboard: vi.fn() }
})

describe('App resilience and accessibility', () => {
  beforeEach(() => vi.mocked(getDashboard).mockReset())

  it('explains snapshot fallback and can retry the live API', async () => {
    vi.mocked(getDashboard)
      .mockRejectedValueOnce(new Error('Network unavailable.'))
      .mockResolvedValueOnce(fallbackDashboard)
    render(<App/>)

    expect(await screen.findByRole('alert')).toHaveTextContent('showing the bundled snapshot')
    expect(screen.getByRole('alert')).toHaveTextContent('Network unavailable.')
    await userEvent.click(screen.getByRole('button', { name: 'Retry connection' }))

    await waitFor(() => expect(screen.getByRole('status')).toHaveTextContent('Live API'))
    expect(screen.queryByRole('alert')).not.toBeInTheDocument()
  })

  it('shows a useful empty search state and clears it', async () => {
    vi.mocked(getDashboard).mockResolvedValue(fallbackDashboard)
    render(<App/>)
    const search = screen.getByRole('textbox', { name: 'Search candidates' })
    await userEvent.type(search, 'not-a-real-compound')

    expect(screen.getByText(/No candidates match/)).toBeInTheDocument()
    await userEvent.click(screen.getByRole('button', { name: 'Clear search' }))
    expect(search).toHaveValue('')
    expect(screen.queryByText(/No candidates match/)).not.toBeInTheDocument()
  })

  it('closes evidence details with Escape and restores page scrolling', async () => {
    vi.mocked(getDashboard).mockResolvedValue(fallbackDashboard)
    render(<App/>)
    await userEvent.click(screen.getByRole('button', { name: `Open ${fallbackDashboard.candidates[0].name}` }))

    expect(screen.getByRole('dialog')).toBeInTheDocument()
    expect(document.body.style.overflow).toBe('hidden')
    fireEvent.keyDown(window, { key: 'Escape' })

    await waitFor(() => expect(screen.queryByRole('dialog')).not.toBeInTheDocument())
    expect(document.body.style.overflow).toBe('')
  })
})
