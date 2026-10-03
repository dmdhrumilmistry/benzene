// SVG renderer for Benzene documents. Produces SVG markup strings in
// "document pixel" space: model coordinates multiplied by style.bondLength.

import { elementColor } from '../core/elements.js';
import { analyzeRings } from '../core/rings.js';
import { add, sub, mul, norm, perp, len, dot, cross, angleOf } from './geometry.js';

export const DEFAULT_STYLE = Object.freeze({
  bondLength: 40,
  lineWidth: 1.6,
  fontSize: 14,
  fontFamily: 'Arial, Helvetica, sans-serif',
  doubleGap: 0.18, // fraction of bond length
  colorAtoms: true,
  showCarbons: false,
  showAtomNumbers: false,
  showValenceErrors: true,
  bondColor: '#1a1a1a',
  selectionColor: '#3b82f6',
  hoverColor: '#f59e0b',
  errorColor: '#e11d48',
  background: null,
});

// Approximate Arial advance widths (em units) for label layout.
const CHAR_W = {
  A: 0.67, B: 0.67, C: 0.72, D: 0.72, E: 0.67, F: 0.61, G: 0.78, H: 0.72, I: 0.28, J: 0.5, K: 0.67, L: 0.56,
  M: 0.83, N: 0.72, O: 0.78, P: 0.67, Q: 0.78, R: 0.72, S: 0.67, T: 0.61, U: 0.72, V: 0.67, W: 0.94, X: 0.67,
  Y: 0.67, Z: 0.61, a: 0.56, b: 0.56, c: 0.5, d: 0.56, e: 0.56, f: 0.28, g: 0.56, h: 0.56, i: 0.22, j: 0.22,
  k: 0.5, l: 0.22, m: 0.83, n: 0.56, o: 0.56, p: 0.56, q: 0.56, r: 0.33, s: 0.5, t: 0.28, u: 0.56, v: 0.5,
  w: 0.72, x: 0.5, y: 0.5, z: 0.5, '+': 0.58, '-': 0.33, '−': 0.58, '•': 0.35, '(': 0.33, ')': 0.33,
};
export function textWidth(s, fs) {
  let w = 0;
  for (const ch of s) w += (CHAR_W[ch] ?? 0.56) * fs;
  return w;
}

export const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const f = (n) => (Math.round(n * 100) / 100).toString();

function chargeText(charge) {
  if (!charge) return '';
  const n = Math.abs(charge);
  return (n > 1 ? n : '') + (charge > 0 ? '+' : '−');
}

/** Direction (unit vector) pointing away from an atom's neighbours. */
export function freeDirection(mol, atomId) {
  const a = mol.getAtom(atomId);
  const nbs = mol.neighbors(atomId);
  if (!nbs.length) return { x: 0.7071, y: -0.7071 };
  let s = { x: 0, y: 0 };
  for (const nb of nbs) {
    const b = mol.getAtom(nb);
    s = add(s, norm(sub(b, a)));
  }
  if (len(s) < 1e-3) {
    // Symmetric (e.g. linear): use the perpendicular of the first bond.
    const b = mol.getAtom(nbs[0]);
    return norm(perp(sub(b, a)));
  }
  return mul(norm(s), -1);
}

/**
 * Decide how an atom is labelled. Returns null for a plain skeletal vertex,
 * otherwise { main, hCount, hSide: 'left'|'right'|'none', charge, isotope, radical, vertexCharge }.
 */
export function atomLabelInfo(mol, atom, style) {
  const degree = mol.degree(atom.id);
  const h = mol.implicitH(atom.id);
  const isC = atom.el === 'C';
  if (atom.label) return { main: atom.label, custom: true, hCount: 0, hSide: 'none', charge: atom.charge, isotope: null, radical: atom.radical };
  if (isC && !style.showCarbons && degree >= 1 && !atom.isotope && !(degree === 1 && atom.charge)) {
    if (atom.charge || atom.radical) return { main: null, vertexCharge: true, charge: atom.charge, radical: atom.radical };
    return null;
  }
  let hSide = 'right';
  if (degree === 0) {
    hSide = ['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I'].includes(atom.el) ? 'left' : 'right';
  } else {
    let dx = 0, dy = 0;
    for (const nb of mol.neighbors(atom.id)) {
      const b = mol.getAtom(nb);
      const v = norm(sub(b, atom));
      dx += v.x; dy += v.y;
    }
    dx /= degree; dy /= degree;
    if (dx > 0.25 || (degree === 1 && dx > 0.05 && Math.abs(dy) < 0.95)) hSide = 'left';
    if (degree >= 2 && Math.abs(dx) < 0.25 && h) hSide = dy > 0.5 ? 'up' : dy < -0.5 ? 'down' : 'right';
  }
  return { main: atom.el, hCount: h, hSide: h ? hSide : 'none', charge: atom.charge, isotope: atom.isotope, radical: atom.radical };
}

/**
 * Render a document to SVG markup (inner content, no <svg> wrapper).
 * opts: { selection: {atoms:Set,bonds:Set,objects:Set}, hover: {type,id}, highlight: {atoms:Set,bonds:Set,color},
 *         annotations: Map<atomId,string>, preview: markup string }
 */
export function renderDocument(doc, style = DEFAULT_STYLE, opts = {}) {
  const s = { ...DEFAULT_STYLE, ...style };
  const mol = doc.mol;
  const bl = s.bondLength;
  const fs = s.fontSize;
  const P = (a) => ({ x: a.x * bl, y: a.y * bl });
  const out = { under: [], bonds: [], atoms: [], over: [] };
  const ringInfo = mol.bondCount ? analyzeRings(mol) : null;
  const sel = opts.selection || { atoms: new Set(), bonds: new Set(), objects: new Set() };
  const hl = opts.highlight;

  // ---- label geometry (for clipping bonds)
  const labels = new Map();
  for (const atom of mol.atoms.values()) {
    const info = atomLabelInfo(mol, atom, s);
    if (!info || !info.main) { labels.set(atom.id, info ? { info, hw: 0, hh: 0 } : null); continue; }
    const mw = textWidth(info.main, fs);
    labels.set(atom.id, { info, hw: (info.custom ? Math.min(mw, fs * 0.8) : mw) / 2 + fs * 0.12, hh: fs * 0.5 + fs * 0.08, mw });
  }
  const clipDist = (atomId, u) => {
    const l = labels.get(atomId);
    if (!l || !l.hw) return 0;
    const tx = Math.abs(u.x) > 1e-6 ? l.hw / Math.abs(u.x) : Infinity;
    const ty = Math.abs(u.y) > 1e-6 ? l.hh / Math.abs(u.y) : Infinity;
    return Math.min(tx, ty);
  };
  const hasLabel = (id) => !!labels.get(id)?.hw;

  // ---- highlights / selection beneath structure
  const selColor = s.selectionColor;
  for (const b of mol.bonds.values()) {
    const p = P(mol.getAtom(b.a1)), q = P(mol.getAtom(b.a2));
    if (hl?.bonds?.has(b.id)) out.under.push(`<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}" stroke="${hl.color || '#fde047'}" stroke-width="${f(bl * 0.28)}" stroke-linecap="round" opacity="0.6"/>`);
    if (sel.bonds?.has(b.id)) out.under.push(`<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}" stroke="${selColor}" stroke-width="${f(bl * 0.22)}" stroke-linecap="round" opacity="0.35"/>`);
    if (opts.hover?.type === 'bond' && opts.hover.id === b.id) out.under.push(`<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}" stroke="${s.hoverColor}" stroke-width="${f(bl * 0.25)}" stroke-linecap="round" opacity="0.45"/>`);
  }
  for (const a of mol.atoms.values()) {
    const p = P(a);
    if (hl?.atoms?.has(a.id)) out.under.push(`<circle cx="${f(p.x)}" cy="${f(p.y)}" r="${f(bl * 0.32)}" fill="${hl.color || '#fde047'}" opacity="0.6"/>`);
    if (sel.atoms?.has(a.id)) out.under.push(`<circle cx="${f(p.x)}" cy="${f(p.y)}" r="${f(bl * 0.25)}" fill="${selColor}" opacity="0.3"/>`);
    if (opts.hover?.type === 'atom' && opts.hover.id === a.id) out.under.push(`<circle cx="${f(p.x)}" cy="${f(p.y)}" r="${f(bl * 0.28)}" fill="${s.hoverColor}" opacity="0.45"/>`);
    if (s.showValenceErrors && mol.hasValenceError(a.id)) out.under.push(`<circle cx="${f(p.x)}" cy="${f(p.y)}" r="${f(bl * 0.34)}" fill="none" stroke="${s.errorColor}" stroke-width="1.2" stroke-dasharray="3 2"/>`);
  }

  // ---- bonds
  const stroke = s.bondColor;
  const lw = s.lineWidth;
  const gap = s.doubleGap * bl;
  const line = (p, q, extra = '') => `<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}" stroke="${stroke}" stroke-width="${f(lw)}" stroke-linecap="round"${extra}/>`;

  const ringCenterFor = (bond) => {
    if (!ringInfo) return null;
    const idxs = ringInfo.bondRings.get(bond.id);
    if (!idxs || !idxs.length) return null;
    // Prefer aromatic rings, then smaller rings, then rings with more double bonds.
    const score = (ring) => {
      let dbl = 0;
      for (let k = 0; k < ring.length; k++) {
        const bb = mol.bondBetween(ring[k], ring[(k + 1) % ring.length]);
        if (bb && bb.order === 2) dbl++;
      }
      return (ringInfo.aromatic.aromaticRings.includes(ring) ? 100 : 0) + dbl * 10 - ring.length;
    };
    const ring = idxs.map((i) => ringInfo.rings[i]).sort((r1, r2) => score(r2) - score(r1))[0];
    let c = { x: 0, y: 0 };
    for (const id of ring) c = add(c, mol.getAtom(id));
    return mul(c, 1 / ring.length);
  };

  const doubleSide = (bond, a, b) => {
    if (bond.position === 'center') return 0;
    if (bond.position === 'left') return 1;
    if (bond.position === 'right') return -1;
    const v = sub(b, a);
    const rc = ringCenterFor(bond);
    if (rc) return cross(v, sub(rc, a)) >= 0 ? 1 : -1;
    const others = (id, skip) => mol.neighbors(id).filter((n) => n !== skip);
    const na = others(bond.a1, bond.a2), nb = others(bond.a2, bond.a1);
    if (!na.length && !nb.length) return 0;
    if ((!na.length && hasLabel(bond.a1)) || (!nb.length && hasLabel(bond.a2))) {
      const rest = na.length ? na : nb;
      if (rest.length !== 1) return 0;
    }
    if ((!na.length && nb.length >= 2) || (!nb.length && na.length >= 2)) return 0;
    let tally = 0;
    for (const n of [...na, ...nb]) {
      const c = cross(v, sub(mol.getAtom(n), a));
      tally += c > 0.01 ? 1 : c < -0.01 ? -1 : 0;
    }
    if (tally === 0) return na.length + nb.length > 0 ? 1 : 0;
    return tally > 0 ? 1 : -1;
  };

  for (const bond of mol.bonds.values()) {
    const A = mol.getAtom(bond.a1), B = mol.getAtom(bond.a2);
    let p = P(A), q = P(B);
    const L = len(sub(q, p));
    if (L < 1e-6) continue;
    const u = mul(sub(q, p), 1 / L);
    const n = perp(u);
    const c1 = clipDist(bond.a1, u), c2 = clipDist(bond.a2, u);
    const p0 = add(p, mul(u, c1)), q0 = sub(q, mul(u, c2));
    const parts = [];
    const st = bond.stereo;
    if (bond.order === 1 && st === 'wedge') {
      const w = gap * 0.75;
      const qq = q0;
      parts.push(`<polygon points="${f(p0.x)},${f(p0.y)} ${f(qq.x + n.x * w)},${f(qq.y + n.y * w)} ${f(qq.x - n.x * w)},${f(qq.y - n.y * w)}" fill="${stroke}" stroke="${stroke}" stroke-width="${f(lw * 0.4)}" stroke-linejoin="round"/>`);
    } else if (bond.order === 1 && st === 'hash') {
      const w = gap * 0.75;
      const segLen = len(sub(q0, p0));
      const count = Math.max(3, Math.round(segLen / (lw * 2.6)));
      for (let i = 0; i <= count; i++) {
        const t = i / count;
        const c = add(p0, mul(sub(q0, p0), t));
        const hw = Math.max(lw * 0.5, w * t);
        parts.push(`<line x1="${f(c.x + n.x * hw)}" y1="${f(c.y + n.y * hw)}" x2="${f(c.x - n.x * hw)}" y2="${f(c.y - n.y * hw)}" stroke="${stroke}" stroke-width="${f(lw * 0.9)}"/>`);
      }
    } else if (bond.order === 1 && st === 'wavy') {
      const segLen = len(sub(q0, p0));
      const waves = Math.max(2, Math.round(segLen / (bl * 0.16)));
      const step = segLen / waves;
      let d = `M${f(p0.x)},${f(p0.y)}`;
      for (let i = 1; i <= waves; i++) {
        const e = add(p0, mul(u, step * i));
        d += ` A${f(step / 2)},${f(step / 2)} 0 0 ${i % 2} ${f(e.x)},${f(e.y)}`;
      }
      parts.push(`<path d="${d}" fill="none" stroke="${stroke}" stroke-width="${f(lw)}"/>`);
    } else if (bond.order === 1 && st === 'bold') {
      parts.push(line(p0, q0, ` stroke-width="${f(gap * 0.55)}" stroke-linecap="butt"`).replace(/stroke-width="[^"]+" stroke-linecap="round"/, ''));
    } else if (bond.order === 1 && st === 'dashed') {
      parts.push(line(p0, q0, ` stroke-dasharray="${f(lw * 2.5)} ${f(lw * 2)}"`));
    } else if (bond.order === 1) {
      parts.push(line(p0, q0));
    } else if (bond.order === 2) {
      const side = doubleSide(bond, A, B);
      if (side === 0) {
        for (const sgn of [1, -1]) {
          const o = mul(n, (sgn * gap) / 2);
          const pp = add(p, o), qq = add(q, o);
          parts.push(line(add(pp, mul(u, clipParallel(bond.a1, u, o, labels))), sub(qq, mul(u, clipParallel(bond.a2, u, o, labels)))));
        }
      } else {
        parts.push(line(p0, q0));
        const o = mul(n, side * gap);
        const pp = add(p, o), qq = add(q, o);
        const sh = gap * 0.9;
        const s1 = hasLabel(bond.a1) ? clipParallel(bond.a1, u, o, labels) : (mol.degree(bond.a1) > 1 ? sh : 0);
        const s2 = hasLabel(bond.a2) ? clipParallel(bond.a2, u, o, labels) : (mol.degree(bond.a2) > 1 ? sh : 0);
        parts.push(line(add(pp, mul(u, s1)), sub(qq, mul(u, s2))));
      }
      if (st === 'wavy' || st === 'dashed') parts[parts.length - 1] = parts[parts.length - 1].replace('/>', ` stroke-dasharray="${f(lw * 2.5)} ${f(lw * 2)}"/>`);
    } else if (bond.order === 3) {
      parts.push(line(p0, q0));
      for (const sgn of [1, -1]) {
        const o = mul(n, sgn * gap);
        parts.push(line(add(add(p, o), mul(u, clipParallel(bond.a1, u, o, labels))), sub(add(q, o), mul(u, clipParallel(bond.a2, u, o, labels)))));
      }
    } else {
      parts.push(line(p0, q0, ' stroke-dasharray="4 3"'));
    }
    out.bonds.push(`<g data-bond="${bond.id}">${parts.join('')}</g>`);
  }

  // ---- atom labels
  for (const atom of mol.atoms.values()) {
    const l = labels.get(atom.id);
    const p = P(atom);
    const color = s.colorAtoms ? elementColor(atom.el) : s.bondColor;
    if (l && l.info.vertexCharge) {
      const d = freeDirection(mol, atom.id);
      const c = add(p, mul(d, fs * 0.75));
      let txt = chargeText(atom.charge);
      if (atom.radical) txt = '•' + txt;
      out.atoms.push(`<text x="${f(c.x)}" y="${f(c.y + fs * 0.32)}" font-size="${f(fs * 0.85)}" text-anchor="middle" fill="${s.bondColor}">${esc(txt)}</text>`);
      continue;
    }
    if (!l || !l.hw) continue;
    const info = l.info;
    const y = p.y + fs * 0.36;
    const parts = [];
    const bg = s.background;
    if (bg) parts.push(`<rect x="${f(p.x - l.hw)}" y="${f(p.y - l.hh)}" width="${f(l.hw * 2)}" height="${f(l.hh * 2)}" fill="${bg}"/>`);
    if (info.custom) {
      const first = textWidth(info.main[0], fs);
      parts.push(`<text x="${f(p.x - first / 2)}" y="${f(y)}" fill="${color}">${esc(info.main)}</text>`);
    } else {
      parts.push(`<text x="${f(p.x)}" y="${f(y)}" text-anchor="middle" fill="${color}">${esc(info.main)}</text>`);
      const mw = l.mw;
      let rightEdge = p.x + mw / 2;
      if (info.hCount) {
        const hw = textWidth('H', fs);
        const sub = info.hCount > 1 ? String(info.hCount) : '';
        const subW = textWidth(sub, fs * 0.7);
        if (info.hSide === 'right') {
          const x0 = p.x + mw / 2;
          parts.push(`<text x="${f(x0)}" y="${f(y)}" fill="${color}">H${sub ? `<tspan dy="${f(fs * 0.25)}" font-size="${f(fs * 0.7)}">${sub}</tspan>` : ''}</text>`);
          rightEdge = x0 + hw + subW;
        } else if (info.hSide === 'left') {
          const x0 = p.x - mw / 2 - hw - subW;
          parts.push(`<text x="${f(x0)}" y="${f(y)}" fill="${color}">H${sub ? `<tspan dy="${f(fs * 0.25)}" font-size="${f(fs * 0.7)}">${sub}</tspan>` : ''}</text>`);
        } else {
          const dy = info.hSide === 'up' ? -fs * 0.95 : fs * 0.95;
          parts.push(`<text x="${f(p.x - (sub ? subW / 2 : 0))}" y="${f(y + dy)}" text-anchor="middle" fill="${color}">H${sub ? `<tspan dy="${f(fs * 0.25)}" font-size="${f(fs * 0.7)}">${sub}</tspan>` : ''}</text>`);
        }
      }
      let ctext = chargeText(info.charge);
      if (info.radical) ctext = (info.radical === 2 ? '••' : '•') + ctext;
      if (ctext) parts.push(`<text x="${f(rightEdge + 1)}" y="${f(y - fs * 0.45)}" font-size="${f(fs * 0.75)}" fill="${color}">${esc(ctext)}</text>`);
      if (info.isotope) {
        const iw = textWidth(String(info.isotope), fs * 0.7);
        const leftEdge = info.hSide === 'left' ? p.x - mw / 2 - textWidth('H', fs) - (info.hCount > 1 ? textWidth(String(info.hCount), fs * 0.7) : 0) : p.x - mw / 2;
        parts.push(`<text x="${f(leftEdge - iw - 0.5)}" y="${f(y - fs * 0.45)}" font-size="${f(fs * 0.7)}" fill="${color}">${info.isotope}</text>`);
      }
    }
    out.atoms.push(`<g data-atom="${atom.id}">${parts.join('')}</g>`);
  }

  // ---- annotations (atom numbers, NMR shifts, ...)
  const annotate = (atomId, text, color) => {
    const a = mol.getAtom(atomId);
    if (!a) return;
    const d = freeDirection(mol, atomId);
    const r = hasLabel(atomId) ? fs * 1.15 : fs * 0.8;
    const c = add(P(a), mul(d, r));
    out.over.push(`<text x="${f(c.x)}" y="${f(c.y + fs * 0.25)}" font-size="${f(fs * 0.62)}" text-anchor="middle" fill="${color}">${esc(text)}</text>`);
  };
  if (s.showAtomNumbers) {
    let i = 1;
    for (const a of mol.atoms.values()) annotate(a.id, String(i++), '#64748b');
  }
  if (opts.annotations) for (const [id, text] of opts.annotations) annotate(id, text, '#0f766e');

  // ---- graphical objects
  for (const o of doc.objects) {
    const selected = sel.objects?.has(o.id);
    const hovered = opts.hover?.type === 'object' && opts.hover.id === o.id;
    out.over.push(renderObject(o, s, { selected, hovered }));
  }

  return [
    `<g class="bz-under">${out.under.join('')}</g>`,
    `<g class="bz-bonds">${out.bonds.join('')}</g>`,
    `<g class="bz-atoms" font-family="${esc(s.fontFamily)}" font-size="${f(fs)}">${out.atoms.join('')}</g>`,
    `<g class="bz-over" font-family="${esc(s.fontFamily)}">${out.over.join('')}</g>`,
    opts.preview || '',
  ].join('');
}

function clipParallel(atomId, u, o, labels) {
  // Clip an offset (parallel) line against the label box of an atom.
  const l = labels.get(atomId);
  if (!l || !l.hw) return 0;
  // Point on offset line closest to atom centre is at offset o; solve for the
  // distance along u where the line exits the box.
  const tx = Math.abs(u.x) > 1e-6 ? (l.hw - Math.sign(u.x) * 0 ) / Math.abs(u.x) : Infinity;
  const ty = Math.abs(u.y) > 1e-6 ? l.hh / Math.abs(u.y) : Infinity;
  let t = Math.min(tx, ty);
  // Slight correction based on lateral offset for nicer results.
  const lateral = Math.abs(dot(o, { x: u.y, y: -u.x }));
  t = Math.max(0, t - lateral * 0.15);
  return t;
}

/** Render a graphical object (arrow / text). Coordinates in model units. */
export function renderObject(o, s, { selected = false, hovered = false } = {}) {
  const bl = s.bondLength;
  const color = o.color || s.bondColor;
  const lw = s.lineWidth;
  const halo = selected ? s.selectionColor : hovered ? s.hoverColor : null;
  if (o.type === 'text') {
    const size = s.fontSize * (o.size || 1);
    const x = o.x * bl, y = o.y * bl;
    const lines = String(o.text).split('\n');
    const tsp = lines.map((ln, i) => `<tspan x="${f(x)}" dy="${i ? f(size * 1.2) : 0}">${richText(ln, size)}</tspan>`).join('');
    let haloEl = '';
    if (halo) {
      const w = Math.max(...lines.map((ln) => textWidth(ln, size)));
      haloEl = `<rect x="${f(x - 3)}" y="${f(y - size)}" width="${f(w + 6)}" height="${f(lines.length * size * 1.2 + 4)}" fill="${halo}" opacity="0.18" rx="3"/>`;
    }
    return `<g data-object="${o.id}">${haloEl}<text x="${f(x)}" y="${f(y)}" font-size="${f(size)}" fill="${color}"${o.bold ? ' font-weight="bold"' : ''}>${tsp}</text></g>`;
  }
  if (o.type === 'arrow') {
    const p = { x: o.x1 * bl, y: o.y1 * bl }, q = { x: o.x2 * bl, y: o.y2 * bl };
    const L = len(sub(q, p));
    if (L < 1) return '';
    const u = mul(sub(q, p), 1 / L), n = perp(u);
    const head = bl * 0.32, hw = bl * 0.11;
    const parts = [];
    const ln = (a, b, extra = '') => `<line x1="${f(a.x)}" y1="${f(a.y)}" x2="${f(b.x)}" y2="${f(b.y)}" stroke="${color}" stroke-width="${f(lw)}" stroke-linecap="round"${extra}/>`;
    const fullHead = (tip, dir) => {
      const nn = perp(dir);
      const base = sub(tip, mul(dir, head));
      const back = sub(tip, mul(dir, head * 0.8));
      return `<polygon points="${f(tip.x)},${f(tip.y)} ${f(base.x + nn.x * hw)},${f(base.y + nn.y * hw)} ${f(back.x)},${f(back.y)} ${f(base.x - nn.x * hw)},${f(base.y - nn.y * hw)}" fill="${color}"/>`;
    };
    const halfHead = (tip, dir, side) => {
      const nn = mul(perp(dir), side);
      const base = sub(tip, mul(dir, head));
      return `<polygon points="${f(tip.x)},${f(tip.y)} ${f(base.x + nn.x * hw * 1.6)},${f(base.y + nn.y * hw * 1.6)} ${f(base.x)},${f(base.y)}" fill="${color}"/>`;
    };
    const style = o.style || 'forward';
    if (halo) parts.push(`<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}" stroke="${halo}" stroke-width="${f(bl * 0.3)}" stroke-linecap="round" opacity="0.25"/>`);
    if (style === 'curved' || style === 'fishhook') {
      const k = o.curve ?? 0.4;
      const ctrl = add(add(p, mul(sub(q, p), 0.5)), mul(n, -k * L));
      const dir = norm(sub(q, ctrl));
      const end = sub(q, mul(dir, head * 0.7));
      parts.push(`<path d="M${f(p.x)},${f(p.y)} Q${f(ctrl.x)},${f(ctrl.y)} ${f(end.x)},${f(end.y)}" fill="none" stroke="${color}" stroke-width="${f(lw)}" stroke-linecap="round"/>`);
      if (halo) parts.unshift(`<path d="M${f(p.x)},${f(p.y)} Q${f(ctrl.x)},${f(ctrl.y)} ${f(q.x)},${f(q.y)}" fill="none" stroke="${halo}" stroke-width="${f(bl * 0.3)}" opacity="0.25"/>`);
      parts.push(style === 'curved' ? fullHead(q, dir) : halfHead(q, dir, cross(sub(ctrl, p), sub(q, p)) > 0 ? 1 : -1));
    } else if (style === 'equilibrium') {
      const o1 = mul(n, -bl * 0.07), o2 = mul(n, bl * 0.07);
      parts.push(ln(add(p, o1), sub(add(q, o1), mul(u, head * 0.6))), halfHead(add(q, o1), u, -1));
      parts.push(ln(add(add(p, o2), mul(u, head * 0.6)), add(q, o2)), halfHead(add(p, o2), mul(u, -1), -1));
    } else if (style === 'resonance') {
      parts.push(ln(add(p, mul(u, head * 0.7)), sub(q, mul(u, head * 0.7))), fullHead(q, u), fullHead(p, mul(u, -1)));
    } else if (style === 'retro') {
      const g = bl * 0.08;
      const tipBack = sub(q, mul(u, head * 0.9));
      parts.push(ln(add(p, mul(n, g)), add(tipBack, mul(n, g))), ln(sub(p, mul(n, g)), sub(tipBack, mul(n, g))));
      parts.push(`<polyline points="${f(tipBack.x + n.x * head * 0.6 - u.x * head * 0.4)},${f(tipBack.y + n.y * head * 0.6 - u.y * head * 0.4)} ${f(q.x)},${f(q.y)} ${f(tipBack.x - n.x * head * 0.6 - u.x * head * 0.4)},${f(tipBack.y - n.y * head * 0.6 - u.y * head * 0.4)}" fill="none" stroke="${color}" stroke-width="${f(lw)}" stroke-linejoin="round" stroke-linecap="round"/>`);
    } else if (style === 'noreaction') {
      parts.push(ln(p, sub(q, mul(u, head * 0.6))), fullHead(q, u));
      const m = add(p, mul(sub(q, p), 0.5));
      const a = rotateVec(u, Math.PI / 4), b = rotateVec(u, -Math.PI / 4);
      parts.push(ln(add(m, mul(a, bl * 0.18)), sub(m, mul(a, bl * 0.18))), ln(add(m, mul(b, bl * 0.18)), sub(m, mul(b, bl * 0.18))));
    } else if (style === 'dashed') {
      parts.push(ln(p, sub(q, mul(u, head * 0.6)), ` stroke-dasharray="${f(lw * 4)} ${f(lw * 3)}"`), fullHead(q, u));
    } else if (style === 'line') {
      parts.push(ln(p, q));
    } else {
      parts.push(ln(p, sub(q, mul(u, head * 0.6))), fullHead(q, u));
    }
    return `<g data-object="${o.id}">${parts.join('')}</g>`;
  }
  return '';
}

function rotateVec(v, t) {
  const c = Math.cos(t), s = Math.sin(t);
  return { x: v.x * c - v.y * s, y: v.x * s + v.y * c };
}

/** Text with simple markup: _x or _{xy} subscript, ^x or ^{xy} superscript. */
export function richText(text, size) {
  let outStr = '';
  const re = /([_^])(\{([^}]*)\}|(\S))/g;
  let last = 0;
  let m;
  while ((m = re.exec(text))) {
    outStr += esc(text.slice(last, m.index));
    const content = m[3] ?? m[4];
    const isSub = m[1] === '_';
    const dy = isSub ? size * 0.25 : -size * 0.4;
    outStr += `<tspan dy="${f(dy)}" font-size="${f(size * 0.7)}">${esc(content)}</tspan><tspan dy="${f(-dy)}">&#8203;</tspan>`;
    last = m.index + m[0].length;
  }
  return outStr + esc(text.slice(last));
}

/** Angle helper re-exported for tools that need label-aware geometry. */
export { angleOf };

/**
 * Full standalone SVG for export. Returns { svg, width, height }.
 * `padding` in px. Embeds the document JSON in <metadata> so the exported
 * image can be re-opened as an editable structure.
 */
export function exportSVG(doc, style = DEFAULT_STYLE, { padding = 12, embed = true, background = null, scale = 1 } = {}) {
  const s = { ...DEFAULT_STYLE, ...style };
  const inner = renderDocument(doc, s, {});
  const b = doc.bbox() || { minX: 0, minY: 0, maxX: 0, maxY: 0 };
  const bl = s.bondLength;
  const extra = s.fontSize * 1.2;
  const minX = b.minX * bl - padding - extra, minY = b.minY * bl - padding - extra;
  const width = (b.maxX - b.minX) * bl + 2 * (padding + extra);
  const height = (b.maxY - b.minY) * bl + 2 * (padding + extra);
  const meta = embed ? `<metadata><benzene:document xmlns:benzene="https://github.com/dmdhrumilmistry/benzene"><![CDATA[${JSON.stringify(doc.toJSON())}]]></benzene:document></metadata>` : '';
  const bg = background ? `<rect x="${f(minX)}" y="${f(minY)}" width="${f(width)}" height="${f(height)}" fill="${background}"/>` : '';
  const svg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="${f(width * scale)}" height="${f(height * scale)}" viewBox="${f(minX)} ${f(minY)} ${f(width)} ${f(height)}">${meta}${bg}${inner}</svg>`;
  return { svg, width: width * scale, height: height * scale };
}

/** Extract an embedded Benzene document from SVG text (or null). */
export function extractDocFromSVG(text) {
  const m = /<benzene:document[^>]*><!\[CDATA\[([\s\S]*?)\]\]><\/benzene:document>/.exec(text);
  if (!m) return null;
  try { return JSON.parse(m[1]); } catch { return null; }
}
