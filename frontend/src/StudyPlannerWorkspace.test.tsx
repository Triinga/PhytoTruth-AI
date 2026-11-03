import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { describe, expect, it } from 'vitest'
import { StudyPlannerWorkspace } from './StudyPlannerWorkspace'

describe('StudyPlannerWorkspace', () => {
  it('shows verified uncertainty and disagreement limits', () => {
    render(<StudyPlannerWorkspace/>)
    expect(screen.getAllByText('0.826')).toHaveLength(2)
    expect(screen.getByText(/95% study-cluster CI 0.375–0.958/)).toBeInTheDocument()
    expect(screen.getByText(/Three disagreement events are insufficient/)).toBeInTheDocument()
    expect(screen.getByText(/reduced MAE from/)).toHaveTextContent('0.517')
  })

  it('switches between outcome-blind experiment recommendations', async () => {
    render(<StudyPlannerWorkspace/>)
    await userEvent.click(screen.getByRole('button', { name: /CHEMBL668/ }))
    expect(screen.getByText('Protocol-complete independent IC50 replication', { selector: 'h3' })).toBeInTheDocument()
    expect(screen.getByText('92%')).toBeInTheDocument()
  })
})
