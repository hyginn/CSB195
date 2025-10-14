/* ----- js/ca.js ----- */
/**
 * Elementary Cellular Automaton core logic (no DOM).
 * Mapping convention: neighborhood 000..111 maps to indices 0..7 respectively.
 * Wolfram code bit i corresponds to neighborhood with index i (LSB = 000, MSB = 111).
 */
;(function (global) {
  'use strict'

  /** Convert integer 0..255 into rule bit array [b0..b7] */
  function decodeRule(code) {
    const bits = new Uint8Array(8)
    for (let i = 0; i < 8; i++) bits[i] = (code >> i) & 1
    return bits // index = neighborhood (0..7)
  }

  /** Convert bit array [b0..b7] back to integer (0..255). */
  function encodeRule(bits) {
    let code = 0
    for (let i = 0; i < 8; i++) if (bits[i]) code |= (1 << i)
    return code >>> 0
  }

  /** Normalize a freestyle 0/1 string into an array of 0/1 numbers. */
  function parseBitString(s) {
    if (!s) return []
    const arr = []
    for (const ch of s) {
      if (ch === '0' || ch === '1') arr.push(ch.charCodeAt(0) & 1)
    }
    return arr
  }

  /** Center-pad or trim an array to target width with zeros. */
  function fitToWidth(bits, width) {
    const out = new Uint8Array(width)
    if (!bits || bits.length === 0) {
      out[Math.floor(width / 2)] = 1
      return out
    }
    if (bits.length >= width) {
      const start = Math.floor((bits.length - width) / 2)
      for (let i = 0; i < width; i++) out[i] = bits[start + i] & 1
      return out
    }
    // center
    const offset = Math.floor((width - bits.length) / 2)
    for (let i = 0; i < bits.length; i++) out[offset + i] = bits[i] & 1
    return out
  }

  /** Compute next generation; returns [row, ruleIndexRow]. */
  function nextGeneration(prev, ruleBits, periodic) {
    const n = prev.length
    const row = new Uint8Array(n)
    const trace = new Uint8Array(n) // which neighborhood index fired
    for (let i = 0; i < n; i++) {
      const L = i === 0 ? (periodic ? prev[n - 1] : 0) : prev[i - 1]
      const C = prev[i]
      const R = i === n - 1 ? (periodic ? prev[0] : 0) : prev[i + 1]
      const idx = (L << 2) | (C << 1) | R // 0..7
      const v = ruleBits[idx]
      row[i] = v
      trace[i] = idx
    }
    return [row, trace]
  }

  /** Run full automaton, returning generations matrix (Uint8Array[]). */
  function runAutomaton(initial, ruleBits, height, periodic) {
    const gens = new Array(height)
    gens[0] = initial
    let prev = initial
    for (let g = 1; g < height; g++) {
      const [row] = nextGeneration(prev, ruleBits, periodic)
      gens[g] = row
      prev = row
    }
    return gens
  }

  /** Random row of width n (0/1, p=0.5). */
  function randomRow(n) {
    const r = new Uint8Array(n)
    for (let i = 0; i < n; i++) r[i] = Math.random() < 0.5 ? 0 : 1
    return r
  }

  // Expose API
  global.ECA = {
    decodeRule,
    encodeRule,
    parseBitString,
    fitToWidth,
    nextGeneration,
    runAutomaton,
    randomRow,
  }
})(window)
