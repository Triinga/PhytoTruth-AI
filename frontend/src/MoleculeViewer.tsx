import { PointerEvent as ReactPointerEvent, WheelEvent as ReactWheelEvent, useEffect, useRef, useState } from 'react'
import { createPortal } from 'react-dom'
import { Atom, Focus, Maximize2, Minus, Plus, RotateCcw, X } from 'lucide-react'
import { fetchMoleculeSvg } from './api'

const MIN_SCALE = 0.5
const MAX_SCALE = 32
const ZOOM_STEP = 1.25
const PRESENTATION_SCALE = 1.6

type ViewBox = { x: number; y: number; width: number; height: number }
type SvgAsset = { markup: string; initialViewBox: ViewBox }

function centeredViewBox(initial: ViewBox, scale: number): ViewBox {
  const width = initial.width / scale
  const height = initial.height / scale
  return {
    x: initial.x + (initial.width - width) / 2,
    y: initial.y + (initial.height - height) / 2,
    width,
    height,
  }
}

function parseDimension(value: string | null, fallback: number) {
  const parsed = Number.parseFloat(value ?? '')
  return Number.isFinite(parsed) && parsed > 0 ? parsed : fallback
}

export function prepareMoleculeSvg(markup: string): SvgAsset {
  const documentNode = new DOMParser().parseFromString(markup, 'image/svg+xml')
  if (documentNode.querySelector('parsererror') || documentNode.documentElement.localName !== 'svg') {
    throw new Error('RDKit returned malformed SVG markup.')
  }
  const root = documentNode.documentElement
  root.querySelectorAll('script,foreignObject,iframe,object,embed,style').forEach(node => node.remove())
  root.querySelectorAll('*').forEach(node => {
    for (const attribute of [...node.attributes]) {
      const name = attribute.name.toLowerCase()
      if (name.startsWith('on')) node.removeAttribute(attribute.name)
      if ((name === 'href' || name.endsWith(':href')) && !attribute.value.startsWith('#')) node.removeAttribute(attribute.name)
      if (name === 'style') {
        const safeStyle = attribute.value
          .replace(/url\s*\([^)]*\)/gi, 'none')
          .replace(/#000000\b/gi, '#dceee8')
          .replace(/#000\b/gi, '#dceee8')
        node.setAttribute(attribute.name, safeStyle)
      }
      if ((name === 'fill' || name === 'stroke') && /^(#000(?:000)?|black)$/i.test(attribute.value)) {
        node.setAttribute(attribute.name, '#dceee8')
      }
    }
  })

  const viewBoxParts = (root.getAttribute('viewBox') ?? '').trim().split(/[ ,]+/).map(Number)
  const initialViewBox: ViewBox = viewBoxParts.length === 4 && viewBoxParts.every(Number.isFinite) && viewBoxParts[2] > 0 && viewBoxParts[3] > 0
    ? { x: viewBoxParts[0], y: viewBoxParts[1], width: viewBoxParts[2], height: viewBoxParts[3] }
    : {
      x: 0,
      y: 0,
      width: parseDimension(root.getAttribute('width'), 520),
      height: parseDimension(root.getAttribute('height'), 260),
    }
  root.setAttribute('viewBox', `${initialViewBox.x} ${initialViewBox.y} ${initialViewBox.width} ${initialViewBox.height}`)
  root.setAttribute('width', '100%')
  root.setAttribute('height', '100%')
  root.setAttribute('preserveAspectRatio', 'xMidYMid meet')
  root.setAttribute('focusable', 'false')
  root.setAttribute('aria-hidden', 'true')
  root.classList.add('rdkit-vector-svg')
  return { markup: new XMLSerializer().serializeToString(root), initialViewBox }
}

function VectorCanvas({ asset, expanded, onExpand }: { asset: SvgAsset; expanded?: boolean; onExpand?: () => void }) {
  const defaultScale = expanded ? 1 : PRESENTATION_SCALE
  const [viewBox, setViewBox] = useState(() => centeredViewBox(asset.initialViewBox, defaultScale))
  const host = useRef<HTMLDivElement | null>(null)
  const drag = useRef<{ pointerId: number; clientX: number; clientY: number; viewBox: ViewBox } | null>(null)
  const initial = asset.initialViewBox
  const scale = initial.width / viewBox.width

  function fitView() {
    drag.current = null
    setViewBox(initial)
  }

  useEffect(() => {
    drag.current = null
    setViewBox(centeredViewBox(asset.initialViewBox, defaultScale))
  }, [asset, defaultScale])

  useEffect(() => {
    const svg = host.current?.querySelector('svg')
    svg?.setAttribute('viewBox', `${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}`)
  }, [asset, viewBox])

  function zoomBy(factor: number, clientX?: number, clientY?: number) {
    setViewBox(current => {
      const currentScale = initial.width / current.width
      const nextScale = Math.min(MAX_SCALE, Math.max(MIN_SCALE, currentScale * factor))
      const nextWidth = initial.width / nextScale
      const nextHeight = initial.height / nextScale
      const bounds = host.current?.getBoundingClientRect()
      const ratioX = clientX !== undefined && bounds?.width ? (clientX - bounds.left) / bounds.width : 0.5
      const ratioY = clientY !== undefined && bounds?.height ? (clientY - bounds.top) / bounds.height : 0.5
      return {
        x: current.x + (current.width - nextWidth) * ratioX,
        y: current.y + (current.height - nextHeight) * ratioY,
        width: nextWidth,
        height: nextHeight,
      }
    })
  }

  function handleWheel(event: ReactWheelEvent<HTMLDivElement>) {
    event.preventDefault()
    event.stopPropagation()
    zoomBy(event.deltaY < 0 ? ZOOM_STEP : 1 / ZOOM_STEP, event.clientX, event.clientY)
  }

  function handlePointerDown(event: ReactPointerEvent<HTMLDivElement>) {
    if (event.button !== 0) return
    event.preventDefault()
    event.stopPropagation()
    event.currentTarget.setPointerCapture?.(event.pointerId)
    drag.current = { pointerId: event.pointerId, clientX: event.clientX, clientY: event.clientY, viewBox }
  }

  function handlePointerMove(event: ReactPointerEvent<HTMLDivElement>) {
    const start = drag.current
    const bounds = host.current?.getBoundingClientRect()
    if (!start || start.pointerId !== event.pointerId || !bounds?.width || !bounds.height) return
    event.preventDefault()
    setViewBox({
      ...start.viewBox,
      x: start.viewBox.x - (event.clientX - start.clientX) * start.viewBox.width / bounds.width,
      y: start.viewBox.y - (event.clientY - start.clientY) * start.viewBox.height / bounds.height,
    })
  }

  function stopDragging(event: ReactPointerEvent<HTMLDivElement>) {
    if (drag.current?.pointerId === event.pointerId) drag.current = null
  }

  return <div className={`molecule-canvas ${expanded ? 'expanded' : 'prediction-structure'}`}>
    <div
      ref={host}
      className="molecule-svg-host"
      data-testid={expanded ? 'expanded-molecule-viewport' : 'molecule-viewport'}
      onWheel={handleWheel}
      onPointerDown={handlePointerDown}
      onPointerMove={handlePointerMove}
      onPointerUp={stopDragging}
      onPointerCancel={stopDragging}
      aria-label="Vector molecule drawing; drag to pan and use the mouse wheel to zoom"
      dangerouslySetInnerHTML={{ __html: asset.markup }}
    />
    <div className="molecule-controls" role="group" aria-label="Molecule view controls">
      {!expanded && <button type="button" onClick={onExpand} aria-label="Open expanded molecule viewer" title="Expand viewer"><Maximize2/></button>}
      <button type="button" onClick={fitView} aria-label="Fit molecule to view" title="Fit to view"><Focus/></button>
      {expanded && <button type="button" onClick={fitView} aria-label="Reset molecule view" title="Reset view"><RotateCcw/></button>}
      <button type="button" onClick={() => zoomBy(ZOOM_STEP)} disabled={scale >= MAX_SCALE} aria-label="Zoom in molecule"><Plus/></button>
      <button type="button" onClick={() => zoomBy(1 / ZOOM_STEP)} disabled={scale <= MIN_SCALE} aria-label="Zoom out molecule"><Minus/></button>
    </div>
    {!expanded && <span className="molecule-canvas-label">2D STRUCTURE · DRAG TO PAN</span>}
    <span className="molecule-zoom-level" aria-live="polite">{Math.round(scale * 100)}%</span>
  </div>
}

function ExpandedViewer({ asset, smiles, onClose }: { asset: SvgAsset; smiles: string; onClose: () => void }) {
  useEffect(() => {
    const previousOverflow = document.body.style.overflow
    document.body.style.overflow = 'hidden'
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') onClose()
    }
    window.addEventListener('keydown', handleKeyDown)
    return () => {
      window.removeEventListener('keydown', handleKeyDown)
      document.body.style.overflow = previousOverflow
    }
  }, [onClose])

  return createPortal(<div className="molecule-modal-backdrop" onMouseDown={event => {
    if (event.target === event.currentTarget) onClose()
  }}>
    <section className="molecule-modal" role="dialog" aria-modal="true" aria-label={`Expanded molecular structure for ${smiles}`}>
      <button type="button" className="molecule-modal-close" onClick={onClose} aria-label="Close expanded molecule viewer"><X/></button>
      <VectorCanvas asset={asset} expanded/>
    </section>
  </div>, document.body)
}

export function MoleculeViewer({ smiles }: { smiles: string }) {
  const [asset, setAsset] = useState<SvgAsset | null>(null)
  const [failed, setFailed] = useState(false)
  const [expanded, setExpanded] = useState(false)

  useEffect(() => {
    const controller = new AbortController()
    setAsset(null)
    setFailed(false)
    setExpanded(false)
    fetchMoleculeSvg(smiles, controller.signal)
      .then(markup => setAsset(prepareMoleculeSvg(markup)))
      .catch(cause => {
        if (!(cause instanceof DOMException && cause.name === 'AbortError')) setFailed(true)
      })
    return () => controller.abort()
  }, [smiles])

  if (failed) return <div className="prediction-structure"><div className="prediction-structure-fallback"><Atom/><span>RDKit vector rendering unavailable</span></div></div>
  if (!asset) return <div className="prediction-structure"><div className="molecule-vector-loading"><Atom/><span>Loading vector structure…</span></div></div>
  return <>
    <VectorCanvas asset={asset} onExpand={() => setExpanded(true)}/>
    {expanded && <ExpandedViewer asset={asset} smiles={smiles} onClose={() => setExpanded(false)}/>} 
  </>
}
