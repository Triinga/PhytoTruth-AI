import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { describe, expect, it } from 'vitest'
import { ExperimentWorkspace } from './ExperimentWorkspace'

describe('ExperimentWorkspace', () => {
  it('reports the locked multi-target negative result and power boundary', () => {
    render(<ExperimentWorkspace/>)
    expect(screen.getByText('Superiority not demonstrated')).toBeInTheDocument()
    expect(screen.getByText(/does not support a superiority claim/)).toBeInTheDocument()
    expect(screen.getByText('1,235')).toBeInTheDocument()
    expect(screen.getByText(/smallest attainable two-sided exact p-value is 0.007813/)).toBeInTheDocument()
  })

  it('interactively switches target-level metrics', async () => {
    render(<ExperimentWorkspace/>)
    await userEvent.selectOptions(screen.getByLabelText('Metric focus'), 'kendall_tau')
    expect(screen.getByText('-0.009')).toBeInTheDocument()
    expect(screen.getByText(/exact p = 0.062500/)).toBeInTheDocument()
    expect(screen.getByText(/Target wins\/ties\/losses: 1\/0\/7/)).toBeInTheDocument()
  })
})
