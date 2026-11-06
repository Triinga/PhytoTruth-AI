import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { describe, expect, it } from 'vitest'
import { ScientificValidationWorkspace } from './ScientificValidationWorkspace'
import results from './data/ethno-convergence-results.json'
import sensitivity from './data/ethno-convergence-sensitivity.json'
import summary from './data/ethno-plants.json'

describe('ScientificValidationWorkspace', () => {
  it('reports the prespecified verdict from the executed artifact', () => {
    render(<ScientificValidationWorkspace/>)
    expect(results.support_rule.verdict).toBe('NOT SUPPORTED')
    expect(screen.getAllByText('NOT SUPPORTED').length).toBeGreaterThan(0)
    expect(screen.getByText('0 of 3 prespecified conditions met')).toBeInTheDocument()
  })

  it('marks every unmet support condition as not met', () => {
    render(<ScientificValidationWorkspace/>)
    const unmet = Object.values(results.support_rule.conditions).filter(met => !met).length
    expect(unmet).toBe(3)
    expect(screen.getAllByText('not met')).toHaveLength(unmet)
  })

  it('links the protocol commit so prespecification can be checked externally', () => {
    render(<ScientificValidationWorkspace/>)
    const protocolLink = screen.getByRole('link', { name: /482bf0d/ })
    expect(protocolLink).toHaveAttribute(
      'href',
      'https://github.com/Triinga/PhytoTruth-AI/commit/482bf0ded20eeb1d2bdc8841d9bbcae3db43b9b4',
    )
    expect(protocolLink).toHaveAttribute('target', '_blank')
  })

  it('states the out-of-domain share rather than hiding it', () => {
    render(<ScientificValidationWorkspace/>)
    const outside = summary.totals.compounds_scored - summary.totals.compounds_in_domain
    const share = ((100 * outside) / summary.totals.compounds_scored).toFixed(1)
    expect(share).toBe('70.6')
    expect(screen.getByText(new RegExp(`${share}% of this chemistry is outside`))).toBeInTheDocument()
  })

  it('shows the in-domain sensitivity analysis without reframing it as the primary outcome', () => {
    render(<ScientificValidationWorkspace/>)
    expect(sensitivity.convergence_supported).toBe(false)
    expect(screen.getByText(/Restricting to in-domain chemistry does not rescue it/)).toBeInTheDocument()
    expect(screen.getByText(/Convergence remains unsupported/)).toBeInTheDocument()
  })

  it('gives a per-plant readout that separates in-domain from out-of-domain compounds', () => {
    render(<ScientificValidationWorkspace/>)
    const plant = summary.plants.find(entry => entry.species === 'coix lacryma-jobi')!
    expect(screen.getByRole('heading', { name: 'Coix Lacryma-Jobi' })).toBeInTheDocument()
    expect(screen.getByText(`${plant.inDomain} in-domain`)).toBeInTheDocument()
    expect(screen.getByText(`${plant.compounds - plant.inDomain} out-of-domain`)).toBeInTheDocument()
    expect(screen.getByText(/Nothing here is clinical proof/)).toBeInTheDocument()
  })

  it('retains Duke attribution for a plant with a recorded traditional use', () => {
    render(<ScientificValidationWorkspace/>)
    const plant = summary.plants.find(entry => entry.species === 'coix lacryma-jobi')!
    expect(plant.isCase).toBe(true)
    expect(screen.getByText('Recorded traditional use')).toBeInTheDocument()
    expect(screen.getByText(plant.activities!)).toBeInTheDocument()
  })

  it('searches the analysed cohort and switches the selected plant', async () => {
    render(<ScientificValidationWorkspace/>)
    await userEvent.type(screen.getByLabelText(/Search .* analysed plants/), 'eucalyptus globulus')
    await userEvent.click(screen.getByRole('button', { name: /Eucalyptus Globulus/ }))
    expect(screen.getByRole('heading', { name: 'Eucalyptus Globulus' })).toBeInTheDocument()
  })

  it('reports no match rather than an empty list for an unknown plant', async () => {
    render(<ScientificValidationWorkspace/>)
    await userEvent.type(screen.getByLabelText(/Search .* analysed plants/), 'zzzzz')
    expect(screen.getByText('No analysed plant matches that name.')).toBeInTheDocument()
  })

  it('declares the high-potential rule as a model-output label, not measured activity', () => {
    render(<ScientificValidationWorkspace/>)
    expect(summary.high_potential_definition.requires_inside_applicability_domain).toBe(true)
    expect(screen.getByText(/Not measured activity/)).toBeInTheDocument()
  })

  it('lists every limitation declared in the protocol', () => {
    render(<ScientificValidationWorkspace/>)
    expect(results.declared_limitations).toHaveLength(7)
    for (const limitation of results.declared_limitations) {
      expect(screen.getByText(limitation)).toBeInTheDocument()
    }
  })

  it('records an input checksum for every source file', () => {
    render(<ScientificValidationWorkspace/>)
    for (const [name, digest] of Object.entries(results.inputs_sha256)) {
      expect(screen.getByText(name)).toBeInTheDocument()
      expect(screen.getByText(`${digest.slice(0, 16)}…`)).toBeInTheDocument()
    }
  })
})
