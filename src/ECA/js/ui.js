/* ----- js/ui.js ----- */
/** Wire the DOM to ECA logic and renderer. */
;(function () {
  'use strict'

  // --- DOM elements ---
  const ruleInput = document.getElementById('ruleInput')
  const widthInput = document.getElementById('widthInput')
  const heightInput = document.getElementById('heightInput')
  const periodicCheckbox = document.getElementById('periodicCheckbox')
  const initialTextarea = document.getElementById('initialTextarea')
  const ruleGrid = document.getElementById('ruleGrid')
  const canvas = document.getElementById('automatonCanvas')
  const btnRun = document.getElementById('btnRun')
  const btnPause = document.getElementById('btnPause')
  const btnStep = document.getElementById('btnStep')
  const btnRender = document.getElementById('btnRender')
  const btnReset = document.getElementById('btnReset')
  const btnRandomize = document.getElementById('btnRandomize')
  const btnExportJSON = document.getElementById('btnExportJSON')
  const btnExportCSV = document.getElementById('btnExportCSV')

  // --- App state ---
  const state = {
    ruleBits: ECA.decodeRule(Number(ruleInput.value) || 0),
    width: Number(widthInput.value) || 101,
    height: Number(heightInput.value) || 150,
    periodic: periodicCheckbox.checked,
    cellSize: 4,
    colors: Array.from({ length: 8 }, (_, i) => defaultColorForIndex(i)),
    generations: null,
    animationHandle: null,
    currentRow: 0,
  }

  function defaultColorForIndex(i) {
    // bs: built with scale via https://leonardocolor.io/scales.html#
    const palette = ["#363c8c", "#3b5f9e", "#427dae", "#4f99ba", "#62b2c3", "#7ac9cc", "#94ded3", "#b1f2dc"]
    return palette[i]
  }

  // --- UI builders ---
  function updateRuleOutputCell(el, bit) {
    el.textContent = String(bit)
    el.classList.toggle('live', bit === 1)
  }

  // Build 8-wide rule+color grid (columns from 111 → 000)
  function buildRuleGrid() {
    ruleGrid.innerHTML = ''
    for (let p = 7; p >= 0; p--) {
      const col = document.createElement('div')
      col.className = 'rulecol'

      // Pattern label
      const pattern = document.createElement('div')
      pattern.className = 'pattern-big'
      pattern.textContent = p.toString(2).padStart(3, '0')

      // Output (click to toggle)
      const out = document.createElement('div')
      out.className = 'out-big'
      out.dataset.index = String(p)
      updateRuleOutputCell(out, state.ruleBits[p])
      out.addEventListener('click', () => {
        state.ruleBits[p] = state.ruleBits[p] ? 0 : 1
        updateRuleOutputCell(out, state.ruleBits[p])
        syncRuleInputFromBits()
        // Optional: live recompute to see effect immediately
        recomputeAll()
      })

      // Hex input + swatch
      const hex = document.createElement('input')
      hex.type = 'text'
      hex.className = 'hex-small'
      hex.value = state.colors[p]
      hex.placeholder = '#RRGGBB'

      const swatch = document.createElement('div')
      swatch.className = 'swatch'
      swatch.style.background = ECARender.normalizeHex(state.colors[p], '#000000')

      hex.addEventListener('input', () => {
        const v = hex.value.trim().toUpperCase()
        state.colors[p] = v
        swatch.style.background = ECARender.normalizeHex(v, '#000000')
        // Optional: live recolor without recomputing
        if (state.generations) {
          ECARender.renderAutomaton(canvas, state.generations, state.ruleBits, state.colors,
            { cellSize: state.cellSize, periodic: state.periodic })
        }
      })

      col.appendChild(pattern)
      col.appendChild(out)
      col.appendChild(hex)
      col.appendChild(swatch)
      ruleGrid.appendChild(col)
    }
  }

  function syncRuleInputFromBits() {
    ruleInput.value = String(ECA.encodeRule(state.ruleBits))
  }

  function syncBitsFromRuleInput() {
    const val = clampInt(Number(ruleInput.value), 0, 255)
    ruleInput.value = String(val)
    state.ruleBits = ECA.decodeRule(val)
    // repaint rule outputs in the grid
    ruleGrid.querySelectorAll('.out-big').forEach((el) => {
      const idx = Number(el.dataset.index)
      updateRuleOutputCell(el, state.ruleBits[idx])
    })
  }

  // --- Helpers ---
  function clampInt(v, lo, hi) {
    v = Math.round(isFinite(v) ? v : lo)
    return Math.max(lo, Math.min(hi, v))
  }

  function computeInitialRow() {
    const bits = ECA.parseBitString(initialTextarea.value)
    const width = state.width
    return ECA.fitToWidth(bits, width)
  }

  function recomputeAll() {
    state.width = clampInt(Number(widthInput.value), 3, 4096)
    widthInput.value = String(state.width)
    state.height = clampInt(Number(heightInput.value), 1, 4096)
    heightInput.value = String(state.height)
    state.periodic = periodicCheckbox.checked

    const initial = computeInitialRow()
    const gens = ECA.runAutomaton(initial, state.ruleBits, state.height, state.periodic)
    state.generations = gens
    state.currentRow = gens.length - 1

    ECARender.renderAutomaton(
      canvas,
      gens,
      state.ruleBits,
      state.colors,
      { cellSize: state.cellSize, periodic: state.periodic }
    )
  }

  // --- Animation ---
  function runAnimated() {
    stopAnimation()
    state.width = clampInt(Number(widthInput.value), 3, 4096)
    state.height = clampInt(Number(heightInput.value), 1, 4096)
    state.periodic = periodicCheckbox.checked

    const initial = computeInitialRow()
    const ctx = canvas.getContext('2d')
    canvas.width = state.width * state.cellSize
    canvas.height = state.height * state.cellSize

    // first row
    ECARender.drawGeneration(ctx, 0, initial, new Uint8Array(state.width), state.colors, state.cellSize)

    let prev = initial
    let rowIndex = 1

    function step() {
      if (rowIndex >= state.height) { stopAnimation(); return }
      const [row, trace] = ECA.nextGeneration(prev, state.ruleBits, state.periodic)
      ECARender.drawGeneration(ctx, rowIndex, row, trace, state.colors, state.cellSize)
      prev = row
      rowIndex++
      state.animationHandle = requestAnimationFrame(step)
    }
    state.animationHandle = requestAnimationFrame(step)
  }

  function stopAnimation() {
    if (state.animationHandle) cancelAnimationFrame(state.animationHandle)
    state.animationHandle = null
  }

  // --- Exporters ---
  function exportJSON() {
    if (!state.generations) recomputeAll()
    const meta = {
      rule: ECA.encodeRule(state.ruleBits),
      width: state.width,
      height: state.height,
      periodic: state.periodic,
      initial: Array.from(state.generations[0]),
      generations: state.generations.map((g) => Array.from(g)),
    }
    const blob = new Blob([JSON.stringify(meta)], { type: 'application/json' })
    downloadBlob(blob, `eca_rule${meta.rule}_w${meta.width}_h${meta.height}.json`)
  }

  function exportCSV() {
    if (!state.generations) recomputeAll()
    const lines = ['generation,cell_index,state']
    for (let g = 0; g < state.generations.length; g++) {
      const row = state.generations[g]
      for (let i = 0; i < row.length; i++) lines.push(`${g},${i},${row[i]}`)
    }
    const blob = new Blob([lines.join('\n')], { type: 'text/csv' })
    const rule = ECA.encodeRule(state.ruleBits)
    downloadBlob(blob, `eca_rule${rule}_w${state.width}_h${state.height}.csv`)
  }

  function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = filename
    document.body.appendChild(a)
    a.click()
    a.remove()
    URL.revokeObjectURL(url)
  }

  // --- Events ---
  ruleInput.addEventListener('change', () => { syncBitsFromRuleInput(); recomputeAll() })
  widthInput.addEventListener('change', recomputeAll)
  heightInput.addEventListener('change', recomputeAll)
  periodicCheckbox.addEventListener('change', recomputeAll)

  btnRender.addEventListener('click', recomputeAll)
  btnReset.addEventListener('click', () => { initialTextarea.value = ''; recomputeAll() })
  btnRandomize.addEventListener('click', () => {
    state.width = clampInt(Number(widthInput.value), 3, 4096)
    const row = ECA.randomRow(state.width)
    initialTextarea.value = Array.from(row).join('')
    recomputeAll()
  })

  btnRun.addEventListener('click', runAnimated)
  btnPause.addEventListener('click', stopAnimation)
  btnStep.addEventListener('click', () => {
    stopAnimation()
    if (!state.generations) recomputeAll()
    const last = state.generations[state.generations.length - 1]
    const [row, trace] = ECA.nextGeneration(last, state.ruleBits, state.periodic)
    const ctx = canvas.getContext('2d')
    const nextIndex = state.currentRow + 1
    if (nextIndex >= state.height) return
    ECARender.drawGeneration(ctx, nextIndex, row, trace, state.colors, state.cellSize)
    state.generations.push(row)
    state.currentRow = nextIndex
  })

  btnExportJSON.addEventListener('click', exportJSON)
  btnExportCSV.addEventListener('click', exportCSV)

  // --- Init ---
  buildRuleGrid()
  recomputeAll()
})();

