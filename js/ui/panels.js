// Side panels: molecular properties, NMR prediction and the structure library.

import { computeProperties } from '../core/properties.js';
import { writeSmiles } from '../core/smiles.js';
import { predictH1, predictC13, simulateSpectrum, formatH1Report, formatC13Report } from '../core/nmr.js';
import { exportSVG, esc } from '../editor/renderer.js';
import { ChemDocument } from '../editor/document.js';
import { Library } from './storage.js';
import { smilesToName } from './pubchem.js';
import { download, safeFilename } from './fileio.js';

const fmt = (n, d = 4) => (Number.isFinite(n) ? n.toFixed(d) : '—');

/** SMILES for a set of atoms (or null on failure). */
export function smilesFor(mol, atomIds = null) {
  try {
    const m = atomIds ? mol.subMolecule(atomIds) : mol;
    if (m.isEmpty()) return '';
    return writeSmiles(m);
  } catch (e) {
    console.warn('SMILES generation failed', e);
    return null;
  }
}

function debounce(fn, ms) {
  let t;
  return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
}

// ===================================================================== props
export class PropertiesPanel {
  constructor(el, app) {
    this.el = el;
    this.app = app;
    this.nameCache = new Map();
    this.update = debounce(() => this._render(), 120);
  }

  _render() {
    const editor = this.app.editor;
    const mol = editor.doc.mol;
    if (mol.isEmpty()) {
      this.el.innerHTML = '<div class="placeholder">Draw a structure to see its properties.</div>';
      return;
    }
    const selIds = editor.selectedAtomIds();
    const ids = selIds.size ? [...selIds] : null;
    let p;
    try { p = computeProperties(mol, ids); } catch (e) { console.error(e); }
    if (!p) { this.el.innerHTML = '<div class="placeholder">Could not compute properties.</div>'; return; }
    const smiles = smilesFor(mol, ids);
    const errors = [...mol.atoms.keys()].filter((id) => mol.hasValenceError(id));
    const lipinski = p.lipinskiViolations === 0 ? '<span class="badge ok">passes</span>' : `<span class="badge warn">${p.lipinskiViolations} violation${p.lipinskiViolations > 1 ? 's' : ''}</span>`;
    const cached = smiles ? this.nameCache.get(smiles) : null;
    this.el.innerHTML = `
      <div class="section">
        <h4>${ids ? 'Selection' : 'Structure'} ${p.fragments > 1 ? `<span class="badge">${p.fragments} fragments</span>` : ''}</h4>
        <div class="formula-big">${p.formulaHTML}</div>
        <div class="row small" id="prop-name">${cached ? this._nameHTML(cached) : `<button class="btn small" id="btn-name">Look up name (PubChem)</button>`}</div>
      </div>
      ${errors.length ? `<div class="section"><h4>Structure check</h4><span class="badge warn">${errors.length} atom${errors.length > 1 ? 's' : ''} with unusual valence</span> <span class="muted small">(circled in red)</span></div>` : ''}
      <div class="section">
        <h4>Mass</h4>
        <dl class="kv">
          <dt>Molecular weight</dt><dd>${fmt(p.molecularWeight, 3)} g/mol</dd>
          <dt>Exact mass</dt><dd>${fmt(p.exactMass, 4)}</dd>
          <dt>m/z</dt><dd>${fmt(p.mz, 4)}${p.charge ? ` (z=${p.charge})` : ''}</dd>
        </dl>
      </div>
      <div class="section">
        <h4>Elemental analysis</h4>
        <div class="mono small">${p.elementalAnalysis.map((e) => `${e.el}, ${e.percent.toFixed(2)}%`).join('; ')}</div>
      </div>
      <div class="section">
        <h4>Isotope pattern (MS)</h4>
        ${this._isotopeSVG(p.isotopePattern)}
      </div>
      <div class="section">
        <h4>Descriptors</h4>
        <dl class="kv">
          <dt>Heavy atoms</dt><dd>${p.heavyAtoms}</dd>
          <dt>Rings (aromatic)</dt><dd>${p.rings} (${p.aromaticRings})</dd>
          <dt>Degree of unsaturation</dt><dd>${p.degreeOfUnsaturation}</dd>
          <dt>H-bond donors</dt><dd>${p.hBondDonors}</dd>
          <dt>H-bond acceptors</dt><dd>${p.hBondAcceptors}</dd>
          <dt>Rotatable bonds</dt><dd>${p.rotatableBonds}</dd>
          <dt>TPSA</dt><dd>${p.tpsa.toFixed(2)} Å²</dd>
          <dt>Lipinski (MW/HBD/HBA)</dt><dd>${lipinski}</dd>
        </dl>
      </div>
      <div class="section">
        <h4>SMILES</h4>
        <div class="copy-row"><input readonly value="${esc(smiles ?? '(unavailable)')}" id="prop-smiles"><button class="btn small" data-copy="prop-smiles">Copy</button></div>
      </div>`;
    this.el.querySelectorAll('[data-copy]').forEach((b) => b.addEventListener('click', () => {
      const v = this.el.querySelector(`#${b.dataset.copy}`).value;
      navigator.clipboard?.writeText(v).then(() => this.app.toast('Copied'));
    }));
    const btn = this.el.querySelector('#btn-name');
    if (btn && smiles) btn.addEventListener('click', () => this._lookupName(smiles));
  }

  _nameHTML(r) {
    if (r === 'none') return '<span class="muted">Not found on PubChem</span>';
    const parts = [];
    if (r.title) parts.push(`<strong>${esc(r.title)}</strong>`);
    if (r.iupac && r.iupac.toLowerCase() !== (r.title || '').toLowerCase()) parts.push(`<span class="muted">${esc(r.iupac)}</span>`);
    parts.push(`<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${r.cid}" target="_blank" rel="noopener">CID ${r.cid}</a>`);
    return `<div>${parts.join('<br>')}</div>`;
  }

  async _lookupName(smiles) {
    const host = this.el.querySelector('#prop-name');
    host.innerHTML = '<span class="muted">Looking up…</span>';
    try {
      const r = await smilesToName(smiles);
      this.nameCache.set(smiles, r || 'none');
      host.innerHTML = this._nameHTML(r || 'none');
      if (r?.title && this.app.docName() === 'Untitled') this.app.setDocName(r.title);
    } catch (e) {
      host.innerHTML = `<span class="muted">Lookup failed: ${esc(e.message)}</span>`;
    }
  }

  _isotopeSVG(peaks) {
    if (!peaks.length) return '';
    const W = 300, H = 90, padL = 6, padR = 6, padB = 18;
    const lo = Math.floor(peaks[0].mass) - 1, hi = Math.ceil(peaks[peaks.length - 1].mass) + 1;
    const x = (m) => padL + ((m - lo) / (hi - lo)) * (W - padL - padR);
    const y = (a) => H - padB - (a / 100) * (H - padB - 12);
    let s = `<svg viewBox="0 0 ${W} ${H}" class="spectrum" style="height:${H}px"><line x1="${padL}" y1="${H - padB}" x2="${W - padR}" y2="${H - padB}" stroke="#94a3b8"/>`;
    // Group into nominal masses for labelling.
    const shown = peaks.filter((p) => p.abundance >= 1);
    for (const p of peaks) s += `<line x1="${x(p.mass)}" y1="${H - padB}" x2="${x(p.mass)}" y2="${y(p.abundance)}" stroke="#2563eb" stroke-width="2"/>`;
    const labelled = new Set();
    for (const p of shown) {
      const nominal = Math.round(p.mass);
      if (labelled.has(nominal)) continue;
      labelled.add(nominal);
      s += `<text x="${x(p.mass)}" y="${H - 5}" font-size="9" text-anchor="middle" fill="#475569">${p.mass.toFixed(p.mass > 1000 ? 1 : 2)}</text>`;
      s += `<text x="${x(p.mass)}" y="${y(p.abundance) - 3}" font-size="8" text-anchor="middle" fill="#64748b">${p.abundance.toFixed(0)}</text>`;
    }
    return s + '</svg>';
  }
}

// ======================================================================= NMR
export class NmrPanel {
  constructor(el, app) {
    this.el = el;
    this.app = app;
    this.nucleus = '1H';
    this.frequency = 400;
    this.showOnStructure = false;
    this.result = null;
    this.update = debounce(() => this._render(), 200);
    this.active = false;
  }

  setActive(on) {
    this.active = on;
    if (on) this._render();
    else this._clearHighlights();
  }

  _clearHighlights() {
    const ed = this.app.editor;
    ed.highlight = null;
    if (ed.annotations) { ed.annotations = null; }
    ed.render();
  }

  _render() {
    if (!this.active) return;
    const editor = this.app.editor;
    const mol = editor.doc.mol;
    if (mol.isEmpty()) {
      this.el.innerHTML = '<div class="placeholder">Draw a structure to predict its NMR spectrum.</div>';
      this._clearHighlights();
      return;
    }
    let res;
    try {
      res = this.nucleus === '1H' ? predictH1(mol, { frequency: this.frequency }) : predictC13(mol);
    } catch (e) {
      console.error(e);
      this.el.innerHTML = `<div class="placeholder">Prediction failed: ${esc(e.message)}</div>`;
      return;
    }
    this.result = res;
    const freq = this.nucleus === '1H' ? this.frequency : Math.round(this.frequency / 4);
    const report = this.nucleus === '1H'
      ? formatH1Report(res, { frequency: this.frequency, solvent: 'CDCl3' })
      : formatC13Report(res, { frequency: freq, solvent: 'CDCl3' });
    const rows = res.signals.map((s, i) => this.nucleus === '1H'
      ? `<tr data-i="${i}"><td class="num">${s.shift.toFixed(2)}</td><td>${s.multiplicity}${s.J?.length ? ` <span class="muted small">J=${s.J.map((j) => j.toFixed(1)).join(', ')}</span>` : ''}</td><td class="num">${s.nH}H</td><td class="small muted">${esc(s.environment || '')}</td></tr>`
      : `<tr data-i="${i}"><td class="num">${s.shift.toFixed(1)}</td><td>${s.multiplicity}</td><td class="num">${s.count > 1 ? `${s.count}C` : ''}</td><td class="small muted">${esc(s.environment || '')}</td></tr>`).join('');
    this.el.innerHTML = `
      <div class="section row">
        <div class="seg" id="nmr-nucleus"><button data-n="1H" class="${this.nucleus === '1H' ? 'active' : ''}"><sup>1</sup>H</button><button data-n="13C" class="${this.nucleus === '13C' ? 'active' : ''}"><sup>13</sup>C</button></div>
        <select id="nmr-freq" title="Spectrometer frequency (1H)">${[60, 90, 300, 400, 500, 600, 800].map((f) => `<option value="${f}"${f === this.frequency ? ' selected' : ''}>${this.nucleus === '1H' ? f : Math.round(f / 4)} MHz</option>`).join('')}</select>
        <label class="row small"><input type="checkbox" id="nmr-annot"${this.showOnStructure ? ' checked' : ''}> Shifts on structure</label>
      </div>
      <div class="section">${this._spectrumSVG(res)}</div>
      ${res.warnings?.length ? `<div class="section small muted">${res.warnings.map(esc).join('<br>')}</div>` : ''}
      <div class="section">
        <table class="signal-table"><thead><tr><th>δ (ppm)</th><th>${this.nucleus === '1H' ? 'Mult.' : 'Type'}</th><th></th><th>Assignment</th></tr></thead><tbody>${rows}</tbody></table>
      </div>
      <div class="section">
        <h4>Report <button class="btn small" id="nmr-copy">Copy</button></h4>
        <div class="report" id="nmr-report">${esc(report).replace(/(1H|13C) NMR/, (m) => (m.startsWith('1H') ? '<sup>1</sup>H NMR' : '<sup>13</sup>C NMR'))}</div>
      </div>
      <div class="section row">
        <button class="btn small" id="nmr-svg">Download spectrum (SVG)</button>
        <button class="btn small" id="nmr-csv">Peak list (CSV)</button>
      </div>
      <p class="muted small">Estimated with additive increment rules (CDCl<sub>3</sub>). Typical accuracy ±0.3 ppm (<sup>1</sup>H) and ±5 ppm (<sup>13</sup>C).</p>`;

    this.el.querySelectorAll('#nmr-nucleus button').forEach((b) => b.addEventListener('click', () => { this.nucleus = b.dataset.n; this._render(); }));
    this.el.querySelector('#nmr-freq').addEventListener('change', (e) => { this.frequency = Number(e.target.value); this._render(); });
    this.el.querySelector('#nmr-annot').addEventListener('change', (e) => { this.showOnStructure = e.target.checked; this._applyAnnotations(); });
    this.el.querySelector('#nmr-copy').addEventListener('click', () => navigator.clipboard?.writeText(report).then(() => this.app.toast('Report copied')));
    this.el.querySelector('#nmr-svg').addEventListener('click', () => {
      const svg = this.el.querySelector('svg.spectrum').outerHTML.replace('<svg ', '<svg xmlns="http://www.w3.org/2000/svg" ');
      download(safeFilename(`${this.app.docName()}_${this.nucleus}_NMR`, 'svg'), svg, 'image/svg+xml');
    });
    this.el.querySelector('#nmr-csv').addEventListener('click', () => {
      const head = this.nucleus === '1H' ? 'shift_ppm,multiplicity,J_Hz,integration,assignment' : 'shift_ppm,type,count,assignment';
      const lines = res.signals.map((s) => (this.nucleus === '1H'
        ? [s.shift.toFixed(2), s.multiplicity, (s.J || []).map((j) => j.toFixed(1)).join(' '), s.nH, `"${s.environment || ''}"`]
        : [s.shift.toFixed(1), s.multiplicity, s.count, `"${s.environment || ''}"`]).join(','));
      download(safeFilename(`${this.app.docName()}_${this.nucleus}_peaks`, 'csv'), [head, ...lines].join('\n'), 'text/csv');
    });
    this.el.querySelectorAll('tr[data-i]').forEach((tr) => {
      tr.addEventListener('mouseenter', () => this._highlightSignal(Number(tr.dataset.i)));
      tr.addEventListener('mouseleave', () => this._highlightSignal(null));
    });
    const spec = this.el.querySelector('svg.spectrum');
    spec.querySelectorAll('[data-sig]').forEach((g) => {
      g.addEventListener('mouseenter', () => this._highlightSignal(Number(g.dataset.sig)));
      g.addEventListener('mouseleave', () => this._highlightSignal(null));
    });
    this._applyAnnotations();
  }

  _highlightSignal(i) {
    const ed = this.app.editor;
    this.el.querySelectorAll('tr[data-i]').forEach((tr) => tr.classList.toggle('active', Number(tr.dataset.i) === i));
    ed.highlight = i === null ? null : { atoms: new Set(this.result.signals[i].atomIds), color: '#fde047' };
    ed.render();
  }

  _applyAnnotations() {
    const ed = this.app.editor;
    if (this.showOnStructure && this.result) {
      const m = new Map();
      for (const s of this.result.signals) for (const id of s.atomIds) m.set(id, this.nucleus === '1H' ? s.shift.toFixed(2) : s.shift.toFixed(1));
      ed.annotations = m;
    } else ed.annotations = null;
    ed.render();
  }

  _spectrumSVG(res) {
    const W = 316, H = 190, padL = 8, padR = 8, padT = 22, padB = 26;
    const is1H = this.nucleus === '1H';
    const shifts = res.signals.map((s) => s.shift);
    let max = is1H ? Math.max(10, ...shifts.map((s) => s + 0.5)) : Math.max(220, ...shifts.map((s) => s + 10));
    let min = is1H ? Math.min(0, ...shifts.map((s) => s - 0.5)) : Math.min(0, ...shifts.map((s) => s - 10));
    if (is1H) { max = Math.ceil(max); min = Math.floor(min); }
    const x = (ppm) => padL + ((max - ppm) / (max - min)) * (W - padL - padR);
    const base = H - padB;
    let s = `<svg viewBox="0 0 ${W} ${H}" class="spectrum" font-family="Arial, sans-serif">`;
    s += `<rect width="${W}" height="${H}" fill="#fff"/>`;
    // Axis
    s += `<line x1="${padL}" y1="${base}" x2="${W - padR}" y2="${base}" stroke="#334155"/>`;
    const step = is1H ? (max - min > 12 ? 2 : 1) : 20;
    for (let t = Math.ceil(min / step) * step; t <= max; t += step) {
      s += `<line x1="${x(t)}" y1="${base}" x2="${x(t)}" y2="${base + 4}" stroke="#334155"/><text x="${x(t)}" y="${base + 14}" font-size="9" text-anchor="middle" fill="#334155">${t}</text>`;
    }
    s += `<text x="${W - padR}" y="${H - 2}" font-size="9" text-anchor="end" fill="#64748b">δ (ppm)</text>`;
    if (is1H) {
      const peaks = [];
      res.signals.forEach((sig, i) => (sig.peaks?.length ? sig.peaks : [{ shift: sig.shift, intensity: sig.nH }]).forEach((p) => peaks.push({ ...p, sig: i })));
      const spec = simulateSpectrum(peaks, { min, max, points: 1600, lineWidth: 0.006 * (400 / this.frequency) + 0.002 });
      let d = '';
      for (let i = 0; i < spec.x.length; i++) d += `${i ? 'L' : 'M'}${x(spec.x[i]).toFixed(1)},${(base - spec.y[i] * (base - padT)).toFixed(1)}`;
      s += `<path d="${d}" fill="none" stroke="#1d4ed8" stroke-width="1"/>`;
      // Invisible hover targets + integration labels.
      let lastLabelX = -Infinity;
      res.signals.forEach((sig, i) => {
        const cx = x(sig.shift);
        const showLabel = cx - lastLabelX > 14;
        if (showLabel) lastLabelX = cx;
        s += `<g data-sig="${i}" style="cursor:pointer"><rect x="${cx - 6}" y="${padT - 14}" width="12" height="${base - padT + 14}" fill="transparent"/>${showLabel ? `<text x="${cx}" y="${padT - 6}" font-size="8" text-anchor="middle" fill="#0f766e">${sig.nH}H</text>` : ''}</g>`;
      });
    } else {
      const maxCount = Math.max(...res.signals.map((g) => g.count || 1));
      let lastLabelX = -Infinity;
      res.signals.forEach((sig, i) => {
        const cx = x(sig.shift);
        const h = (0.45 + 0.55 * ((sig.count || 1) / maxCount)) * (base - padT);
        const color = sig.multiplicity === 'C' ? '#7c3aed' : '#1d4ed8';
        s += `<g data-sig="${i}" style="cursor:pointer"><rect x="${cx - 4}" y="${padT - 14}" width="8" height="${base - padT + 14}" fill="transparent"/><line x1="${cx}" y1="${base}" x2="${cx}" y2="${base - h}" stroke="${color}" stroke-width="1.4"/>${Math.abs(cx - lastLabelX) > 16 ? `<text x="${cx}" y="${base - h - 3}" font-size="7.5" text-anchor="middle" fill="#475569">${sig.shift.toFixed(1)}</text>` : ''}</g>`;
        if (Math.abs(cx - lastLabelX) > 16) lastLabelX = cx;
      });
    }
    return s + '</svg>';
  }
}

// =================================================================== library
export class LibraryPanel {
  constructor(el, app) {
    this.el = el;
    this.app = app;
    this.filter = '';
  }

  render() {
    const entries = Library.list().filter((e) => !this.filter || e.name.toLowerCase().includes(this.filter) || (e.formula || '').toLowerCase().includes(this.filter));
    this.el.innerHTML = `
      <div class="lib-save"><input type="text" id="lib-name" placeholder="Name" value="${esc(this.app.docName())}"><button class="btn primary" id="lib-save">Save</button></div>
      <div class="row" style="margin-bottom:10px">
        <input type="search" id="lib-filter" placeholder="Filter…" value="${esc(this.filter)}" style="flex:1;min-width:0">
        <button class="btn small" id="lib-export" title="Download library as JSON">Export</button>
        <button class="btn small" id="lib-import" title="Import a library JSON file">Import</button>
      </div>
      ${entries.length ? `<div class="lib-list">${entries.map((e) => this._card(e)).join('')}</div>` : '<div class="placeholder">No saved structures yet.<br>Press <kbd>Ctrl</kbd>+<kbd>S</kbd> to save the current drawing.</div>'}
      <p class="muted small">Saved in this browser only (localStorage). Use Export to back up or move your library.</p>`;
    this.el.querySelector('#lib-save').addEventListener('click', () => this.app.saveToLibrary(this.el.querySelector('#lib-name').value));
    this.el.querySelector('#lib-name').addEventListener('keydown', (e) => { if (e.key === 'Enter') this.app.saveToLibrary(e.target.value); e.stopPropagation(); });
    const filter = this.el.querySelector('#lib-filter');
    filter.addEventListener('input', (e) => { this.filter = e.target.value.toLowerCase(); this.render(); this.el.querySelector('#lib-filter').focus(); });
    filter.addEventListener('keydown', (e) => e.stopPropagation());
    this.el.querySelector('#lib-export').addEventListener('click', () => download('benzene-library.json', Library.exportJSON(), 'application/json'));
    this.el.querySelector('#lib-import').addEventListener('click', () => document.getElementById('library-input').click());
    this.el.querySelectorAll('.lib-card').forEach((card) => {
      const id = card.dataset.id;
      card.querySelector('.lib-thumb').addEventListener('click', () => this.app.openLibraryEntry(id));
      card.querySelector('.lib-meta').addEventListener('click', () => this.app.openLibraryEntry(id));
      card.querySelectorAll('[data-act]').forEach((b) => b.addEventListener('click', (ev) => {
        ev.stopPropagation();
        const act = b.dataset.act;
        if (act === 'insert') this.app.insertLibraryEntry(id);
        else if (act === 'rename') {
          const e = Library.get(id);
          const name = prompt('Rename structure', e.name);
          if (name) { Library.update(id, { name }); this.render(); }
        } else if (act === 'download') {
          const e = Library.get(id);
          download(safeFilename(e.name, 'benzene.json'), JSON.stringify(e.doc, null, 1), 'application/json');
        } else if (act === 'delete') {
          const e = Library.get(id);
          if (confirm(`Delete “${e.name}” from the library?`)) { Library.remove(id); this.render(); }
        }
      }));
    });
  }

  _card(e) {
    const date = new Date(e.updated || e.created).toLocaleDateString();
    return `<div class="lib-card" data-id="${e.id}">
      <div class="lib-thumb" title="Open">${e.thumbnail || ''}</div>
      <div class="lib-meta"><div class="lib-name" title="${esc(e.name)}">${esc(e.name)}</div><div class="lib-sub"><span>${esc(e.formula || '')}</span><span>${date}</span></div></div>
      <div class="lib-actions">
        <button data-act="insert" title="Insert into current drawing">Insert</button>
        <button data-act="rename" title="Rename">Rename</button>
        <button data-act="download" title="Download">⤓</button>
        <button data-act="delete" title="Delete">✕</button>
      </div>
    </div>`;
  }
}

/** Small standalone SVG thumbnail for a document (no metadata). */
export function thumbnailSVG(doc, style) {
  if (doc.isEmpty()) return '';
  const { svg } = exportSVG(doc, { ...style, bondLength: 22, fontSize: 10, lineWidth: 1.1 }, { embed: false, padding: 4 });
  return svg.replace(/^<\?xml[^>]*>\s*/, '');
}

export { ChemDocument };
