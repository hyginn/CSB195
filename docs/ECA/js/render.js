/* ----- js/render.js ----- */
/** Canvas rendering utilities for ECA. */
;(function (global) {
  'use strict'

  /** Validate a hex color; return uppercase normalized string or default. */
  function normalizeHex(hex, fallback) {
    if (typeof hex !== 'string') return fallback
    const s = hex.trim().toUpperCase()
    const ok = /^#([0-9A-F]{6})$/.test(s)
    return ok ? s : fallback
  }

  /** Create ImageData for a single generation row. */
  function drawGeneration(ctx, y, row, trace, colors, cellSize) {
    const n = row.length
    const w = n * cellSize
    // We draw using fillRect for clarity. ImageData is faster for very large grids.
    for (let x = 0; x < n; x++) {
      const alive = row[x]
      if (!alive) {
        ctx.fillStyle = '#FFFFFF'
      } else {
        const idx = trace ? trace[x] : 0
        const color = colors[idx] || '#000000'
        ctx.fillStyle = normalizeHex(color, '#000000')
      }
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize)
    }
  }

  /** Render full automaton given generations and (optional) per-step traces. */
  function renderAutomaton(canvas, generations, ruleBits, colors, options) {
    const cellSize = (options && options.cellSize) || 4
    const height = generations.length
    const width = generations[0].length
    const ctx = canvas.getContext('2d', { willReadFrequently: false })
    canvas.width = width * cellSize
    canvas.height = height * cellSize

    // For coloring by rule, we need to recompute trace row-by-row on the fly.
    // This maintains zero dependency from computation stage if caller didn't keep traces.
    const periodic = options && !!options.periodic
    let prev = generations[0]

    // First row: color by current cell (no rule fired yet) — use #FFFFFF for zeros, and fallback color for ones.
    const firstTrace = new Uint8Array(width) // zeros
    drawGeneration(ctx, 0, prev, firstTrace, colors, cellSize)

    for (let g = 1; g < height; g++) {
      const [row, trace] = global.ECA.nextGeneration(prev, ruleBits, periodic)
      drawGeneration(ctx, g, row, trace, colors, cellSize)
      prev = row
    }
  }

  global.ECARender = { normalizeHex, drawGeneration, renderAutomaton }
})(window)

