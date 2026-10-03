// Interactive structure editor: tools, hit testing, history and rendering.

import { ChemDocument, normalizeMolecule } from './document.js';
import { renderDocument, DEFAULT_STYLE, renderObject, textWidth, freeDirection } from './renderer.js';
import {
  add, sub, mul, dist, norm, perp, angleOf, fromAngle, snapAngle, distToSegment, rotateAround, rectFromPoints,
  pointInRect, cross, rad, mid,
} from './geometry.js';
import { isElement } from '../core/elements.js';
import { ABBREVIATIONS } from './templates.js';

const ATOM_HIT = 0.32;
const BOND_HIT = 0.2;
const MERGE_DIST = 0.22;
const DRAG_THRESHOLD = 0.12;

/** Hotkeys that change the element of the hovered atom. */
const ELEMENT_KEYS = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', h: 'H', i: 'I', l: 'Cl', b: 'Br', k: 'K', m: 'Mg', a: 'Al', g: 'Si', e: 'Se', B: 'B', L: 'Li', N: 'Na' };

export class Editor extends EventTarget {
  constructor(svg, { parseSmiles = null, layoutMolecule = null } = {}) {
    super();
    this.svg = svg;
    this.world = svg.querySelector('#world');
    this.doc = new ChemDocument();
    this.style = { ...DEFAULT_STYLE };
    this.view = { zoom: 1, x: 0, y: 0 };
    this.tool = { name: 'bond', order: 1, stereo: 'none' };
    this.element = 'C';
    this.selection = { atoms: new Set(), bonds: new Set(), objects: new Set() };
    this.hover = null;
    this.highlight = null;
    this.annotations = null;
    this.undoStack = [];
    this.redoStack = [];
    this.clipboard = null;
    this.drag = null;
    this.preview = '';
    this.pointer = { x: 0, y: 0 };
    this.parseSmiles = parseSmiles;
    this.layoutMolecule = layoutMolecule;
    this.spaceDown = false;
    this._bindEvents();
    this.resize();
    this.view.x = this.svg.clientWidth / 2;
    this.view.y = this.svg.clientHeight / 2;
    this.render();
  }

  // ================================================================ basics
  get bl() { return this.style.bondLength; }

  setDocument(doc, { keepHistory = false } = {}) {
    if (keepHistory) this._snapshot();
    else { this.undoStack = []; this.redoStack = []; }
    this.doc = doc;
    this.clearSelection();
    this.fitToView();
    this._changed();
  }

  setTool(tool) {
    this.tool = { ...tool };
    this.preview = '';
    this.drag = null;
    this.dispatchEvent(new CustomEvent('tool', { detail: this.tool }));
    this.render();
  }

  setStyle(patch) {
    Object.assign(this.style, patch);
    this.render();
  }

  toModel(clientX, clientY) {
    const r = this.svg.getBoundingClientRect();
    return {
      x: (clientX - r.left - this.view.x) / this.view.zoom / this.bl,
      y: (clientY - r.top - this.view.y) / this.view.zoom / this.bl,
    };
  }

  toScreen(p) {
    return { x: p.x * this.bl * this.view.zoom + this.view.x, y: p.y * this.bl * this.view.zoom + this.view.y };
  }

  resize() {
    const w = this.svg.clientWidth, h = this.svg.clientHeight;
    this.svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  }

  // ============================================================== history
  _serialize() { return JSON.stringify(this.doc.toJSON()); }

  _snapshot() {
    this.undoStack.push(this._serialize());
    if (this.undoStack.length > 200) this.undoStack.shift();
    this.redoStack = [];
  }

  /** Run a mutation with undo support. */
  change(fn) {
    const before = this._serialize();
    const result = fn();
    if (this._serialize() !== before) {
      this.undoStack.push(before);
      if (this.undoStack.length > 200) this.undoStack.shift();
      this.redoStack = [];
      this._changed();
    } else this.render();
    return result;
  }

  undo() {
    if (!this.undoStack.length) return;
    this.redoStack.push(this._serialize());
    this.doc = ChemDocument.fromJSON(JSON.parse(this.undoStack.pop()));
    this.clearSelection(false);
    this._changed();
  }

  redo() {
    if (!this.redoStack.length) return;
    this.undoStack.push(this._serialize());
    this.doc = ChemDocument.fromJSON(JSON.parse(this.redoStack.pop()));
    this.clearSelection(false);
    this._changed();
  }

  _changed() {
    // Drop selection entries that no longer exist.
    for (const id of this.selection.atoms) if (!this.doc.mol.getAtom(id)) this.selection.atoms.delete(id);
    for (const id of this.selection.bonds) if (!this.doc.mol.getBond(id)) this.selection.bonds.delete(id);
    for (const id of this.selection.objects) if (!this.doc.getObject(id)) this.selection.objects.delete(id);
    if (this.hover && !this._hoverExists()) this.hover = null;
    this.render();
    this.dispatchEvent(new CustomEvent('change'));
  }

  _hoverExists() {
    const h = this.hover;
    if (h.type === 'atom') return !!this.doc.mol.getAtom(h.id);
    if (h.type === 'bond') return !!this.doc.mol.getBond(h.id);
    if (h.type === 'object') return !!this.doc.getObject(h.id);
    return true;
  }

  // ============================================================ selection
  clearSelection(emit = true) {
    this.selection = { atoms: new Set(), bonds: new Set(), objects: new Set() };
    if (emit) { this.render(); this.dispatchEvent(new CustomEvent('selection')); }
  }

  hasSelection() {
    const s = this.selection;
    return s.atoms.size + s.bonds.size + s.objects.size > 0;
  }

  selectAll() {
    this.selection = {
      atoms: new Set(this.doc.mol.atoms.keys()),
      bonds: new Set(this.doc.mol.bonds.keys()),
      objects: new Set(this.doc.objects.map((o) => o.id)),
    };
    this.setTool({ name: 'select' });
    this.dispatchEvent(new CustomEvent('selection'));
  }

  /** Atom ids covered by the selection (selected atoms + atoms of selected bonds). */
  selectedAtomIds() {
    const ids = new Set(this.selection.atoms);
    for (const bid of this.selection.bonds) {
      const b = this.doc.mol.getBond(bid);
      if (b) { ids.add(b.a1); ids.add(b.a2); }
    }
    return ids;
  }

  _selectBondsBetweenSelectedAtoms() {
    for (const b of this.doc.mol.bonds.values()) {
      if (this.selection.atoms.has(b.a1) && this.selection.atoms.has(b.a2)) this.selection.bonds.add(b.id);
    }
  }

  selectionBBox() {
    const ids = this.selectedAtomIds();
    let box = ids.size ? this.doc.mol.bbox(ids) : null;
    for (const oid of this.selection.objects) {
      const o = this.doc.getObject(oid);
      if (!o) continue;
      const pts = o.type === 'arrow' ? [{ x: o.x1, y: o.y1 }, { x: o.x2, y: o.y2 }] : [{ x: o.x, y: o.y - 0.4 }, { x: o.x + this._textWidthModel(o), y: o.y + 0.2 }];
      for (const p of pts) {
        if (!box) box = { minX: p.x, minY: p.y, maxX: p.x, maxY: p.y };
        box.minX = Math.min(box.minX, p.x); box.maxX = Math.max(box.maxX, p.x);
        box.minY = Math.min(box.minY, p.y); box.maxY = Math.max(box.maxY, p.y);
      }
    }
    if (box) { box.cx = (box.minX + box.maxX) / 2; box.cy = (box.minY + box.maxY) / 2; }
    return box;
  }

  _textWidthModel(o) {
    const size = this.style.fontSize * (o.size || 1);
    return Math.max(...String(o.text).split('\n').map((l) => textWidth(l, size))) / this.bl;
  }

  // ========================================================== hit testing
  hitTest(p, { skipAtoms = null } = {}) {
    const mol = this.doc.mol;
    let best = null, bestD = ATOM_HIT;
    for (const a of mol.atoms.values()) {
      if (skipAtoms?.has(a.id)) continue;
      const d = dist(p, a);
      if (d < bestD) { bestD = d; best = { type: 'atom', id: a.id }; }
    }
    if (best) return best;
    bestD = BOND_HIT;
    for (const b of mol.bonds.values()) {
      if (skipAtoms?.has(b.a1) || skipAtoms?.has(b.a2)) continue;
      const d = distToSegment(p, mol.getAtom(b.a1), mol.getAtom(b.a2));
      if (d < bestD) { bestD = d; best = { type: 'bond', id: b.id }; }
    }
    if (best) return best;
    for (let i = this.doc.objects.length - 1; i >= 0; i--) {
      const o = this.doc.objects[i];
      if (o.type === 'arrow') {
        const pts = this._arrowPoints(o);
        for (let k = 0; k < pts.length - 1; k++) if (distToSegment(p, pts[k], pts[k + 1]) < 0.22) return { type: 'object', id: o.id };
      } else if (o.type === 'text') {
        const size = (this.style.fontSize * (o.size || 1)) / this.bl;
        const lines = String(o.text).split('\n').length;
        const r = { minX: o.x - 0.1, minY: o.y - size, maxX: o.x + this._textWidthModel(o) + 0.1, maxY: o.y + (lines - 1) * size * 1.2 + 0.25 * size };
        if (pointInRect(p, r)) return { type: 'object', id: o.id };
      }
    }
    return null;
  }

  _arrowPoints(o) {
    const p = { x: o.x1, y: o.y1 }, q = { x: o.x2, y: o.y2 };
    if (o.style !== 'curved' && o.style !== 'fishhook') return [p, q];
    const L = dist(p, q);
    const n = perp(norm(sub(q, p)));
    const c = add(mid(p, q), mul(n, -(o.curve ?? 0.4) * L));
    const pts = [];
    for (let t = 0; t <= 1.0001; t += 0.1) {
      const a = (1 - t) * (1 - t), b = 2 * (1 - t) * t, cc = t * t;
      pts.push({ x: a * p.x + b * c.x + cc * q.x, y: a * p.y + b * c.y + cc * q.y });
    }
    return pts;
  }

  _atomNear(p, skip = null) {
    let best = null, bd = MERGE_DIST;
    for (const a of this.doc.mol.atoms.values()) {
      if (skip?.has(a.id)) continue;
      const d = dist(p, a);
      if (d < bd) { bd = d; best = a; }
    }
    return best;
  }

  // ===================================================== geometry helpers
  /** Best angle (radians) for a new bond leaving atomId. */
  bestBondAngle(atomId, order = 1) {
    const mol = this.doc.mol;
    const a = mol.getAtom(atomId);
    const nbs = mol.neighbors(atomId);
    if (!nbs.length) return rad(-30);
    const angles = nbs.map((n) => angleOf(sub(mol.getAtom(n), a)));
    let candidates = [];
    if (nbs.length === 1) {
      const t = angles[0];
      const bond = mol.bondsOf(atomId)[0];
      if (bond.order === 3 || order === 3 || (bond.order === 2 && order === 2)) return t + Math.PI;
      const nb = mol.getAtom(nbs[0]);
      const nbOthers = mol.neighbors(nb.id).filter((x) => x !== atomId);
      candidates = [t + rad(120), t - rad(120)];
      if (nbOthers.length) {
        // Zigzag: choose the candidate on the opposite side of the neighbour's other substituent.
        const o = mol.getAtom(nbOthers[0]);
        const v = sub(a, nb);
        const side = Math.sign(cross(v, sub(o, nb)));
        candidates.sort((x, y) => {
          const sx = Math.sign(cross(v, fromAngle(x))), sy = Math.sign(cross(v, fromAngle(y)));
          return (sx === side ? 1 : 0) - (sy === side ? 1 : 0);
        });
      } else {
        // Prefer pointing up-right / right-ish.
        candidates.sort((x, y) => Math.cos(y) - Math.sin(y) * 0.2 - (Math.cos(x) - Math.sin(x) * 0.2));
      }
    } else {
      const sorted = [...angles].sort((x, y) => x - y);
      let bestGap = -1, bestMid = 0;
      for (let i = 0; i < sorted.length; i++) {
        const s0 = sorted[i], s1 = i + 1 < sorted.length ? sorted[i + 1] : sorted[0] + 2 * Math.PI;
        if (s1 - s0 > bestGap) { bestGap = s1 - s0; bestMid = (s0 + s1) / 2; }
      }
      candidates = [bestMid];
    }
    // Avoid collisions with existing atoms.
    for (const c of candidates) {
      const q = add(a, fromAngle(c));
      const clash = [...mol.atoms.values()].some((x) => x.id !== atomId && dist(x, q) < 0.4);
      if (!clash) return c;
    }
    let best = candidates[0], bestScore = -Infinity;
    for (let k = 0; k < 24; k++) {
      const c = rad(k * 15);
      const q = add(a, fromAngle(c));
      let m = Infinity;
      for (const x of mol.atoms.values()) if (x.id !== atomId) m = Math.min(m, dist(x, q));
      if (m > bestScore) { bestScore = m; best = c; }
    }
    return best;
  }

  /** Create (or reuse) an atom at point p. */
  _atomAt(p, el = 'C') {
    const existing = this._atomNear(p);
    if (existing) return existing;
    return this.doc.mol.addAtom({ el, x: p.x, y: p.y });
  }

  _applyBondTool(bond) {
    const t = this.tool;
    const mol = this.doc.mol;
    const stereo = t.stereo || 'none';
    if (stereo !== 'none') {
      if (bond.stereo === stereo && bond.order === 1) { const tmp = bond.a1; bond.a1 = bond.a2; bond.a2 = tmp; mol.invalidate(); }
      else { bond.order = 1; bond.stereo = stereo; }
      return;
    }
    if (bond.stereo !== 'none') { bond.stereo = 'none'; bond.order = t.order; return; }
    if (t.order === 1) bond.order = bond.order === 1 ? 2 : bond.order === 2 ? 3 : 1;
    else if (t.order === 2) {
      if (bond.order === 2) {
        const cycle = ['auto', 'center', 'left', 'right'];
        bond.position = cycle[(cycle.indexOf(bond.position) + 1) % cycle.length];
      } else bond.order = 2;
    } else bond.order = t.order;
    if (bond.order !== 2) bond.position = 'auto';
  }

  /** Add a bond from an atom in the best free direction. Returns the new atom. */
  addBondFromAtom(atomId, order = 1, stereo = 'none', el = 'C') {
    const mol = this.doc.mol;
    const a = mol.getAtom(atomId);
    const q = add(a, fromAngle(this.bestBondAngle(atomId, order)));
    const b = this._atomAt(q, el);
    if (b.id === atomId) return b;
    const existing = mol.bondBetween(atomId, b.id);
    if (existing) return b;
    mol.addBond(atomId, b.id, order, stereo);
    return b;
  }

  // ======================================================== ring geometry
  ringPlacement(hit, p, tmpl) {
    const mol = this.doc.mol;
    const n = tmpl.size;
    const R = 1 / (2 * Math.sin(Math.PI / n));
    const apothem = 1 / (2 * Math.tan(Math.PI / n));
    const delta = (2 * Math.PI) / n;
    let center, theta0, dir = 1, fixed = [];
    if (hit?.type === 'bond') {
      const bond = mol.getBond(hit.id);
      const a = mol.getAtom(bond.a1), b = mol.getAtom(bond.a2);
      const u = norm(sub(b, a)), nn = perp(u);
      let tally = 0;
      for (const id of [bond.a1, bond.a2]) {
        for (const nb of mol.neighbors(id)) {
          if (nb === bond.a1 || nb === bond.a2) continue;
          tally += Math.sign(cross(u, sub(mol.getAtom(nb), a)));
        }
      }
      let side = tally > 0 ? -1 : tally < 0 ? 1 : (cross(u, sub(p, a)) >= 0 ? 1 : -1);
      // perp(u) = (-u.y, u.x); cross(u, nn) = 1 so +1 side == +nn
      center = add(mid(a, b), mul(nn, side * apothem));
      theta0 = angleOf(sub(a, center));
      const tb = angleOf(sub(b, center));
      let d = tb - theta0;
      while (d <= -Math.PI) d += 2 * Math.PI;
      while (d > Math.PI) d -= 2 * Math.PI;
      dir = d > 0 ? 1 : -1;
      fixed = [a.id, b.id];
    } else if (hit?.type === 'atom') {
      const a = mol.getAtom(hit.id);
      const t = mol.degree(a.id) ? this.bestBondAngle(a.id) : rad(-90);
      center = add(a, fromAngle(t, R));
      theta0 = angleOf(sub(a, center));
      fixed = [a.id];
    } else {
      center = p;
      theta0 = -Math.PI / 2 + (n % 2 === 0 ? 0 : 0);
      if (n === 4 || n === 8) theta0 += delta / 2;
    }
    const pts = [];
    for (let k = 0; k < n; k++) pts.push(add(center, fromAngle(theta0 + dir * k * delta, R)));
    const ids = pts.map((pt, k) => {
      if (fixed[k] !== undefined) return fixed[k];
      const near = this._atomNear(pt);
      return near ? near.id : null;
    });
    return { pts, ids, center };
  }

  placeRing(hit, p, tmpl) {
    const mol = this.doc.mol;
    const { pts, ids } = this.ringPlacement(hit, p, tmpl);
    const n = pts.length;
    const atomIds = pts.map((pt, k) => ids[k] ?? mol.addAtom({ x: pt.x, y: pt.y }).id);
    const newBonds = [];
    for (let k = 0; k < n; k++) {
      const a = atomIds[k], b = atomIds[(k + 1) % n];
      let bond = mol.bondBetween(a, b);
      if (!bond) { bond = mol.addBond(a, b, 1); newBonds.push(bond); }
    }
    if (tmpl.aromatic) {
      const hasDouble = (id) => mol.bondsOf(id).some((x) => x.order === 2);
      const tryPattern = (start) => {
        const chosen = [];
        const used = new Set();
        for (let k = start; k < n; k += 2) {
          const a = atomIds[k], b = atomIds[(k + 1) % n];
          const bond = mol.bondBetween(a, b);
          if (!newBonds.includes(bond)) continue;
          if (hasDouble(a) || hasDouble(b) || used.has(a) || used.has(b)) continue;
          chosen.push(bond); used.add(a); used.add(b);
        }
        return chosen;
      };
      const maxDoubles = n === 5 ? 2 : Math.floor(n / 2);
      const options = [tryPattern(0), tryPattern(1)].map((c) => c.slice(0, maxDoubles));
      const best = options.sort((x, y) => y.length - x.length)[0];
      for (const b of best) b.order = 2;
    }
    return atomIds;
  }

  // ======================================================= chain geometry
  chainPoints(start, p, startAtomId = null) {
    const v = sub(p, start);
    const L = Math.hypot(v.x, v.y);
    const step = Math.cos(rad(30));
    const nBonds = Math.max(1, Math.round(L / step));
    const theta = snapAngle(angleOf(v), 30);
    let sign = -1;
    if (startAtomId) {
      const mol = this.doc.mol;
      const nbs = mol.neighbors(startAtomId);
      if (nbs.length === 1) {
        const a = mol.getAtom(startAtomId), nb = mol.getAtom(nbs[0]);
        const c = cross(fromAngle(theta), sub(nb, a));
        sign = c > 0 ? -1 : 1;
      }
    }
    const pts = [{ ...start }];
    let cur = { ...start };
    for (let i = 0; i < nBonds; i++) {
      const t = theta + sign * rad(30) * (i % 2 === 0 ? 1 : -1);
      cur = add(cur, fromAngle(t));
      pts.push(cur);
    }
    return pts;
  }

  // ================================================================ edits
  deleteSelection() {
    if (!this.hasSelection()) return;
    this.change(() => {
      for (const id of this.selection.bonds) this.doc.mol.removeBond(id);
      for (const id of this.selection.atoms) this.doc.mol.removeAtom(id);
      for (const id of this.selection.objects) this.doc.removeObject(id);
      this._removeOrphans();
    });
    this.clearSelection();
  }

  _deleteHit(hit) {
    const mol = this.doc.mol;
    if (hit.type === 'atom') mol.removeAtom(hit.id);
    else if (hit.type === 'bond') { mol.removeBond(hit.id); this._removeOrphans(); }
    else if (hit.type === 'object') this.doc.removeObject(hit.id);
  }

  /** Remove isolated carbon atoms left behind after bond deletion. */
  _removeOrphans() {
    const mol = this.doc.mol;
    for (const a of mol.atomList()) {
      if (a.el === 'C' && mol.degree(a.id) === 0 && !a.charge && !a.isotope && !a.label && a._keep !== true) mol.removeAtom(a.id);
    }
  }

  mergeAtoms(keepId, removeId) {
    const mol = this.doc.mol;
    if (keepId === removeId) return;
    for (const b of [...mol.bondsOf(removeId)]) {
      const other = mol.otherAtom(b, removeId);
      if (other === keepId || mol.bondBetween(keepId, other)) continue;
      const a1 = b.a1 === removeId ? keepId : b.a1;
      const a2 = b.a2 === removeId ? keepId : b.a2;
      mol.addBond(a1, a2, b.order, b.stereo, { position: b.position });
    }
    mol.removeAtom(removeId);
  }

  setAtomElement(atomId, el) {
    const a = this.doc.mol.getAtom(atomId);
    if (!a) return;
    a.el = el; a.label = null; a.hCount = null;
    if (el === 'C' || el === 'H') a.charge = a.charge && el === 'C' ? a.charge : 0;
  }

  changeCharge(atomId, delta) {
    const a = this.doc.mol.getAtom(atomId);
    if (a) a.charge = Math.max(-4, Math.min(4, (a.charge || 0) + delta));
  }

  /**
   * Apply typed text to an atom: element symbol ("N"), element with H/charge
   * ("NH2", "O-", "NH3+"), abbreviation ("OMe", "CO2H") or arbitrary label (R, Ar, X).
   */
  applyAtomLabel(atomId, text) {
    const mol = this.doc.mol;
    const atom = mol.getAtom(atomId);
    text = text.trim();
    if (!atom || !text) return;
    const m = /^(?:H(\d*))?([A-Z][a-z]?)(?:H(\d*))?(\d*[+-]|[+-]\d*)?$/.exec(text);
    if (m && isElement(m[2]) && !(m[2] === 'C' && (m[1] !== undefined || m[3] !== undefined) && ABBREVIATIONS[text])) {
      atom.el = m[2]; atom.label = null; atom.hCount = null;
      const cs = m[4] || '';
      if (cs) {
        const num = parseInt(cs.replace(/[+-]/g, ''), 10) || 1;
        atom.charge = cs.includes('-') ? -num : num;
      } else atom.charge = 0;
      return;
    }
    const abbr = ABBREVIATIONS[text] ?? Object.entries(ABBREVIATIONS).find(([k]) => k.toLowerCase() === text.toLowerCase())?.[1];
    if (abbr && this.parseSmiles) {
      this.expandGroup(atomId, abbr);
      return;
    }
    atom.el = 'R'; atom.label = text; atom.hCount = 0; atom.charge = 0;
  }

  /** Replace an atom with a group given by SMILES (first atom = attachment point). */
  expandGroup(atomId, smiles) {
    const mol = this.doc.mol;
    const atom = mol.getAtom(atomId);
    const g = this.parseSmiles(smiles);
    normalizeMolecule(g);
    const gAtoms = g.atomList();
    const g0 = gAtoms[0];
    // Direction from attachment atom into the group.
    let into = { x: 0, y: 0 };
    for (const a of gAtoms) if (a.id !== g0.id) into = add(into, sub(a, g0));
    if (Math.hypot(into.x, into.y) < 1e-6) into = { x: 1, y: 0 };
    let desired;
    const nbs = mol.neighbors(atomId);
    if (nbs.length) {
      let s = { x: 0, y: 0 };
      for (const n of nbs) s = add(s, norm(sub(atom, mol.getAtom(n))));
      desired = Math.hypot(s.x, s.y) > 1e-6 ? s : { x: 1, y: 0 };
    } else desired = { x: 1, y: 0 };
    const rot = angleOf(desired) - angleOf(into);
    const map = new Map();
    for (const a of gAtoms) {
      const r = rotateAround(a, g0, rot);
      const pos = { x: atom.x + (r.x - g0.x), y: atom.y + (r.y - g0.y) };
      if (a.id === g0.id) {
        atom.el = a.el; atom.charge = a.charge; atom.hCount = a.hCount; atom.isotope = a.isotope; atom.label = null;
        map.set(a.id, atom.id);
      } else {
        const near = this._atomNear(pos);
        map.set(a.id, near && !nbs.includes(near.id) ? near.id : mol.addAtom({ ...a, id: undefined, x: pos.x, y: pos.y }).id);
      }
    }
    for (const b of g.bonds.values()) {
      const a1 = map.get(b.a1), a2 = map.get(b.a2);
      if (a1 !== a2 && !mol.bondBetween(a1, a2)) mol.addBond(a1, a2, b.order, b.stereo);
    }
  }

  /** Insert a molecule (e.g. from SMILES or a file) centred at point p (model units). */
  insertMolecule(m, p = null) {
    const center = p ?? this._viewCenterModel();
    normalizeMolecule(m, center.x, center.y);
    let ids;
    this.change(() => {
      const map = this.doc.mol.merge(m);
      ids = [...map.values()];
    });
    this.selection = { atoms: new Set(ids), bonds: new Set(), objects: new Set() };
    this._selectBondsBetweenSelectedAtoms();
    this.render();
    this.dispatchEvent(new CustomEvent('selection'));
    return ids;
  }

  _viewCenterModel() {
    const r = this.svg.getBoundingClientRect();
    const c = this.toModel(r.left + r.width / 2, r.top + r.height / 2);
    // Avoid overlapping existing content: place to the right of it.
    const b = this.doc.bbox();
    if (b && !this.doc.isEmpty()) {
      const inside = c.x > b.minX - 1 && c.x < b.maxX + 1 && c.y > b.minY - 1 && c.y < b.maxY + 1;
      if (inside) return { x: b.maxX + 3, y: b.cy };
    }
    return c;
  }

  /** Re-layout (clean up) selected fragments, or everything. */
  cleanStructure() {
    if (!this.layoutMolecule) return;
    const mol = this.doc.mol;
    const sel = this.selectedAtomIds();
    this.change(() => {
      for (const frag of mol.fragments()) {
        if (sel.size && !frag.some((id) => sel.has(id))) continue;
        const sub = mol.subMolecule(frag);
        const before = mol.bbox(frag);
        this.layoutMolecule(sub);
        const after = sub.bbox();
        for (const id of frag) {
          const a = sub.getAtom(id), t = mol.getAtom(id);
          t.x = a.x - after.cx + before.cx;
          t.y = a.y - after.cy + before.cy;
        }
      }
    });
  }

  rotateSelection(angle) {
    const box = this.selectionBBox() || this.doc.bbox();
    if (!box) return;
    const c = { x: box.cx, y: box.cy };
    const all = !this.hasSelection();
    this.change(() => this._transformSelection((p) => rotateAround(p, c, angle), all));
  }

  flipSelection(axis) {
    const box = this.selectionBBox() || this.doc.bbox();
    if (!box) return;
    const all = !this.hasSelection();
    this.change(() => {
      this._transformSelection((p) => (axis === 'h' ? { x: 2 * box.cx - p.x, y: p.y } : { x: p.x, y: 2 * box.cy - p.y }), all);
      // Mirroring inverts stereo: swap wedge <-> hash for selected bonds.
      const ids = all ? new Set(this.doc.mol.atoms.keys()) : this.selectedAtomIds();
      for (const b of this.doc.mol.bonds.values()) {
        if (!ids.has(b.a1) || !ids.has(b.a2)) continue;
        if (b.stereo === 'wedge') b.stereo = 'hash';
        else if (b.stereo === 'hash') b.stereo = 'wedge';
      }
    });
  }

  _transformSelection(fn, all = false) {
    const ids = all ? new Set(this.doc.mol.atoms.keys()) : this.selectedAtomIds();
    for (const id of ids) {
      const a = this.doc.mol.getAtom(id);
      const q = fn(a); a.x = q.x; a.y = q.y;
    }
    const objs = all ? this.doc.objects : this.doc.objects.filter((o) => this.selection.objects.has(o.id));
    for (const o of objs) {
      if (o.type === 'arrow') {
        const p1 = fn({ x: o.x1, y: o.y1 }), p2 = fn({ x: o.x2, y: o.y2 });
        o.x1 = p1.x; o.y1 = p1.y; o.x2 = p2.x; o.y2 = p2.y;
      } else {
        const q = fn({ x: o.x, y: o.y }); o.x = q.x; o.y = q.y;
      }
    }
  }

  nudge(dx, dy) {
    if (!this.hasSelection()) return;
    this.change(() => this._transformSelection((p) => ({ x: p.x + dx, y: p.y + dy })));
  }

  // ============================================================ clipboard
  copySelection() {
    if (!this.hasSelection()) return null;
    const ids = this.selectedAtomIds();
    const sub = this.doc.mol.subMolecule(ids);
    const objects = this.doc.objects.filter((o) => this.selection.objects.has(o.id)).map((o) => ({ ...o }));
    this.clipboard = { molecule: sub.toJSON(), objects };
    return { mol: sub, objects };
  }

  cutSelection() {
    const r = this.copySelection();
    this.deleteSelection();
    return r;
  }

  paste(data = this.clipboard) {
    if (!data) return;
    const offset = 0.8;
    const newAtoms = [], newObjs = [];
    this.change(() => {
      const map = new Map();
      for (const a of data.molecule.atoms) {
        const { id, ...rest } = a;
        const na = this.doc.mol.addAtom({ ...rest, x: a.x + offset, y: a.y + offset });
        map.set(id, na.id); newAtoms.push(na.id);
      }
      for (const b of data.molecule.bonds) {
        if (map.has(b.a1) && map.has(b.a2)) this.doc.mol.addBond(map.get(b.a1), map.get(b.a2), b.order, b.stereo, { position: b.position });
      }
      for (const o of data.objects || []) {
        const { id, ...rest } = o;
        const no = { ...rest };
        if (no.type === 'arrow') { no.x1 += offset; no.y1 += offset; no.x2 += offset; no.y2 += offset; }
        else { no.x += offset; no.y += offset; }
        newObjs.push(this.doc.addObject(no).id);
      }
    });
    // Shift clipboard so repeated pastes cascade.
    this.clipboard = {
      molecule: { atoms: data.molecule.atoms.map((a) => ({ ...a, x: a.x + offset, y: a.y + offset })), bonds: data.molecule.bonds },
      objects: (data.objects || []).map((o) => (o.type === 'arrow' ? { ...o, x1: o.x1 + offset, y1: o.y1 + offset, x2: o.x2 + offset, y2: o.y2 + offset } : { ...o, x: o.x + offset, y: o.y + offset })),
    };
    this.selection = { atoms: new Set(newAtoms), bonds: new Set(), objects: new Set(newObjs) };
    this._selectBondsBetweenSelectedAtoms();
    this.setTool({ name: 'select' });
    this.dispatchEvent(new CustomEvent('selection'));
  }

  // ================================================================= view
  zoomBy(factor, cx = null, cy = null) {
    const r = this.svg.getBoundingClientRect();
    const sx = cx ?? r.width / 2, sy = cy ?? r.height / 2;
    const z = Math.max(0.15, Math.min(8, this.view.zoom * factor));
    const k = z / this.view.zoom;
    this.view.x = sx - (sx - this.view.x) * k;
    this.view.y = sy - (sy - this.view.y) * k;
    this.view.zoom = z;
    this.render();
    this.dispatchEvent(new CustomEvent('view'));
  }

  fitToView() {
    const b = this.doc.bbox();
    const w = this.svg.clientWidth, h = this.svg.clientHeight;
    if (!b) { this.view = { zoom: 1, x: w / 2, y: h / 2 }; this.render(); this.dispatchEvent(new CustomEvent('view')); return; }
    const bw = (b.width + 3) * this.bl, bh = (b.height + 3) * this.bl;
    const z = Math.max(0.2, Math.min(2, Math.min(w / bw, h / bh)));
    this.view.zoom = z;
    this.view.x = w / 2 - b.cx * this.bl * z;
    this.view.y = h / 2 - b.cy * this.bl * z;
    this.render();
    this.dispatchEvent(new CustomEvent('view'));
  }

  // =============================================================== render
  render() {
    const t = `translate(${this.view.x},${this.view.y}) scale(${this.view.zoom})`;
    this.world.setAttribute('transform', t);
    const overlay = this._overlayMarkup();
    this.world.innerHTML = renderDocument(this.doc, this.style, {
      selection: this.selection,
      hover: this.drag?.kind === 'move' ? null : this.hover,
      highlight: this.highlight,
      annotations: this.annotations,
      preview: this.preview + overlay,
    });
  }

  _overlayMarkup() {
    const bl = this.bl, z = this.view.zoom;
    let s = '';
    if (this.drag?.kind === 'rubber' && this.drag.current) {
      const r = rectFromPoints(this.drag.start, this.drag.current);
      s += `<rect x="${r.minX * bl}" y="${r.minY * bl}" width="${(r.maxX - r.minX) * bl}" height="${(r.maxY - r.minY) * bl}" fill="rgba(59,130,246,0.08)" stroke="#3b82f6" stroke-width="${1 / z}" stroke-dasharray="${4 / z} ${3 / z}"/>`;
    }
    if (this.tool.name === 'select' && this.hasSelection() && this.drag?.kind !== 'rubber') {
      const b = this.selectionBBox();
      if (b) {
        const pad = 0.35;
        const x = (b.minX - pad) * bl, y = (b.minY - pad) * bl, w = (b.maxX - b.minX + 2 * pad) * bl, h = (b.maxY - b.minY + 2 * pad) * bl;
        const hx = (b.cx) * bl, hy = (b.minY - pad - 0.7) * bl;
        s += `<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="none" stroke="#3b82f6" stroke-width="${1 / z}" stroke-dasharray="${5 / z} ${3 / z}"/>`;
        s += `<line x1="${hx}" y1="${y}" x2="${hx}" y2="${hy}" stroke="#3b82f6" stroke-width="${1 / z}"/>`;
        s += `<circle cx="${hx}" cy="${hy}" r="${6 / z}" fill="#fff" stroke="#3b82f6" stroke-width="${1.5 / z}"/>`;
        this._rotateHandle = { x: b.cx, y: b.minY - pad - 0.7, cx: b.cx, cy: b.cy };
      }
    } else this._rotateHandle = null;
    return s;
  }

  _ghost(pts, closed = false, ids = null) {
    const bl = this.bl;
    const c = '#3b82f6';
    let s = '';
    for (let i = 0; i < pts.length - (closed ? 0 : 1); i++) {
      const p = pts[i], q = pts[(i + 1) % pts.length];
      s += `<line x1="${p.x * bl}" y1="${p.y * bl}" x2="${q.x * bl}" y2="${q.y * bl}" stroke="${c}" stroke-width="2" stroke-linecap="round" opacity="0.55"/>`;
    }
    pts.forEach((p, i) => {
      if (!ids || ids[i] === null) s += `<circle cx="${p.x * bl}" cy="${p.y * bl}" r="2.5" fill="${c}" opacity="0.55"/>`;
    });
    return s;
  }

  // =============================================================== events
  _bindEvents() {
    const svg = this.svg;
    svg.addEventListener('pointerdown', (e) => this._onDown(e));
    svg.addEventListener('pointermove', (e) => this._onMove(e));
    svg.addEventListener('pointerup', (e) => this._onUp(e));
    svg.addEventListener('pointerleave', () => { if (!this.drag) { this.hover = null; this.preview = ''; this.render(); } });
    svg.addEventListener('dblclick', (e) => this._onDblClick(e));
    svg.addEventListener('wheel', (e) => {
      e.preventDefault();
      const r = svg.getBoundingClientRect();
      if (e.ctrlKey || e.metaKey) this.zoomBy(Math.exp(-e.deltaY * 0.01), e.clientX - r.left, e.clientY - r.top);
      else { this.view.x -= e.deltaX; this.view.y -= e.deltaY; this.render(); }
    }, { passive: false });
    svg.addEventListener('contextmenu', (e) => e.preventDefault());
    window.addEventListener('resize', () => { this.resize(); this.render(); });
  }

  _onDown(e) {
    if (this._textEditor) return;
    const p = this.toModel(e.clientX, e.clientY);
    this.svg.setPointerCapture(e.pointerId);
    if (e.button === 1 || e.button === 2 || this.spaceDown || this.tool.name === 'pan') {
      this.drag = { kind: 'pan', sx: e.clientX, sy: e.clientY, vx: this.view.x, vy: this.view.y };
      return;
    }
    const hit = this.hitTest(p);
    const t = this.tool.name;
    const base = { start: p, hit, moved: false, shift: e.shiftKey, alt: e.altKey };

    if (t === 'select') {
      if (this._rotateHandle && dist(p, this._rotateHandle) < 12 / this.view.zoom / this.bl) {
        this.drag = { ...base, kind: 'rotate', center: { x: this._rotateHandle.cx, y: this._rotateHandle.cy }, a0: angleOf(sub(p, { x: this._rotateHandle.cx, y: this._rotateHandle.cy })), snap: this._serialize(), orig: this.doc.toJSON() };
        return;
      }
      if (hit) {
        const key = hit.type === 'atom' ? 'atoms' : hit.type === 'bond' ? 'bonds' : 'objects';
        if (!this.selection[key].has(hit.id)) {
          if (!e.shiftKey) this.selection = { atoms: new Set(), bonds: new Set(), objects: new Set() };
          this.selection[key].add(hit.id);
          if (hit.type === 'bond') { const b = this.doc.mol.getBond(hit.id); this.selection.atoms.add(b.a1); this.selection.atoms.add(b.a2); }
          this.dispatchEvent(new CustomEvent('selection'));
        } else if (e.shiftKey) {
          this.selection[key].delete(hit.id);
          this.dispatchEvent(new CustomEvent('selection'));
        }
        this.drag = { ...base, kind: 'move', last: p, snap: this._serialize() };
      } else {
        if (!e.shiftKey) this.clearSelection();
        this.drag = { ...base, kind: 'rubber', current: p, prevSel: e.shiftKey ? this._copySel() : null };
      }
      this.render();
      return;
    }
    if (t === 'erase') {
      this.drag = { ...base, kind: 'erase', snap: this._serialize(), erased: false };
      if (hit) { this._deleteHit(hit); this.drag.erased = true; this.hover = null; this.render(); }
      return;
    }
    this.drag = { ...base, kind: t };
  }

  _copySel() {
    return { atoms: new Set(this.selection.atoms), bonds: new Set(this.selection.bonds), objects: new Set(this.selection.objects) };
  }

  _onMove(e) {
    const p = this.toModel(e.clientX, e.clientY);
    this.pointer = p;
    this.dispatchEvent(new CustomEvent('pointer', { detail: p }));
    const d = this.drag;
    if (!d) { this._updateHover(p); return; }
    if (d.kind === 'pan') {
      this.view.x = d.vx + (e.clientX - d.sx);
      this.view.y = d.vy + (e.clientY - d.sy);
      this.render();
      return;
    }
    if (!d.moved && dist(p, d.start) > DRAG_THRESHOLD) d.moved = true;
    if (!d.moved) return;
    const mol = this.doc.mol;

    switch (d.kind) {
      case 'move': {
        const dx = p.x - d.last.x, dy = p.y - d.last.y;
        d.last = p;
        this._transformSelection((q) => ({ x: q.x + dx, y: q.y + dy }));
        this.render();
        break;
      }
      case 'rotate': {
        let ang = angleOf(sub(p, d.center)) - d.a0;
        if (!e.altKey) ang = snapAngle(ang, 15);
        this.doc = ChemDocument.fromJSON(d.orig);
        this._transformSelection((q) => rotateAround(q, d.center, ang));
        this.preview = `<text x="${p.x * this.bl + 14}" y="${p.y * this.bl - 10}" font-size="12" fill="#3b82f6" font-family="sans-serif">${Math.round((ang * 180) / Math.PI)}°</text>`;
        this.render();
        break;
      }
      case 'rubber': {
        d.current = p;
        const r = rectFromPoints(d.start, p);
        const sel = d.prevSel ? this._copySel() : { atoms: new Set(), bonds: new Set(), objects: new Set() };
        if (d.prevSel) Object.assign(sel, { atoms: new Set(d.prevSel.atoms), bonds: new Set(d.prevSel.bonds), objects: new Set(d.prevSel.objects) });
        for (const a of mol.atoms.values()) if (pointInRect(a, r)) sel.atoms.add(a.id);
        for (const b of mol.bonds.values()) if (sel.atoms.has(b.a1) && sel.atoms.has(b.a2)) sel.bonds.add(b.id);
        for (const o of this.doc.objects) {
          const pt = o.type === 'arrow' ? mid({ x: o.x1, y: o.y1 }, { x: o.x2, y: o.y2 }) : { x: o.x, y: o.y };
          if (pointInRect(pt, r)) sel.objects.add(o.id);
        }
        this.selection = sel;
        this.render();
        break;
      }
      case 'erase': {
        const hit = this.hitTest(p);
        if (hit) { this._deleteHit(hit); d.erased = true; this.render(); }
        break;
      }
      case 'bond':
      case 'atom': {
        const from = d.hit?.type === 'atom' ? mol.getAtom(d.hit.id) : d.start;
        const end = this._bondDragEnd(from, p, e.altKey, d.hit?.type === 'atom' ? d.hit.id : null);
        d.end = end;
        this.preview = this._ghost([from, end.pos]);
        this.render();
        break;
      }
      case 'chain': {
        const startAtom = d.hit?.type === 'atom' ? d.hit.id : null;
        const from = startAtom ? mol.getAtom(startAtom) : d.start;
        const pts = this.chainPoints(from, p, startAtom);
        d.pts = pts;
        const last = pts[pts.length - 1];
        this.preview = this._ghost(pts) + `<text x="${last.x * this.bl + 10}" y="${last.y * this.bl - 8}" font-size="12" fill="#3b82f6" font-family="sans-serif">${pts.length - 1}</text>`;
        this.render();
        break;
      }
      case 'arrow': {
        let q = p;
        const style = this.tool.style || 'forward';
        if (!e.altKey && style !== 'curved' && style !== 'fishhook') {
          const L = dist(d.start, p);
          q = add(d.start, fromAngle(snapAngle(angleOf(sub(p, d.start)), 15), L));
        }
        d.end = q;
        this.preview = renderObject({ type: 'arrow', style, x1: d.start.x, y1: d.start.y, x2: q.x, y2: q.y, color: '#3b82f6' }, this.style);
        this.render();
        break;
      }
      default:
        break;
    }
  }

  _bondDragEnd(from, p, free, fromId) {
    const target = this.hitTest(p);
    if (target?.type === 'atom' && target.id !== fromId) {
      const a = this.doc.mol.getAtom(target.id);
      return { pos: { x: a.x, y: a.y }, atomId: target.id };
    }
    const v = sub(p, from);
    const t = free ? angleOf(v) : snapAngle(angleOf(v), 15);
    const L = free ? Math.max(0.3, Math.hypot(v.x, v.y)) : 1;
    const pos = add(from, fromAngle(t, L));
    const near = this._atomNear(pos, fromId ? new Set([fromId]) : null);
    return near ? { pos: { x: near.x, y: near.y }, atomId: near.id } : { pos, atomId: null };
  }

  _updateHover(p) {
    const hit = this.hitTest(p);
    const changed = JSON.stringify(hit) !== JSON.stringify(this.hover);
    this.hover = hit;
    let preview = '';
    if (this.tool.name === 'ring' && this.tool.template) {
      const { pts, ids } = this.ringPlacement(hit, p, this.tool.template);
      preview = this._ghost(pts, true, ids);
    }
    if (changed || preview || this.preview) {
      this.preview = preview;
      this.render();
      this.dispatchEvent(new CustomEvent('hover', { detail: hit }));
    }
    this.svg.style.cursor = this._cursorFor(hit);
  }

  _cursorFor(hit) {
    const t = this.tool.name;
    if (this.spaceDown || t === 'pan') return 'grab';
    if (t === 'select') return hit ? 'move' : 'default';
    if (t === 'erase') return hit ? 'pointer' : 'default';
    if (t === 'text') return 'text';
    return 'crosshair';
  }

  _onUp(e) {
    const d = this.drag;
    this.drag = null;
    try { this.svg.releasePointerCapture(e.pointerId); } catch { /* ignore */ }
    if (!d) return;
    const p = this.toModel(e.clientX, e.clientY);
    const mol = this.doc.mol;
    this.preview = '';

    switch (d.kind) {
      case 'pan': this.dispatchEvent(new CustomEvent('view')); break;
      case 'move': {
        if (d.moved) {
          // Merge dragged atoms onto nearby stationary atoms.
          const moving = this.selectedAtomIds();
          for (const id of [...moving]) {
            const a = mol.getAtom(id);
            if (!a) continue;
            const near = this._atomNear(a, moving);
            if (near) { this.mergeAtoms(near.id, id); this.selection.atoms.delete(id); }
          }
          this.undoStack.push(d.snap);
          this.redoStack = [];
          this._changed();
        } else if (!d.shift && d.hit) {
          // Plain click on selected item selects only it.
          const key = d.hit.type === 'atom' ? 'atoms' : d.hit.type === 'bond' ? 'bonds' : 'objects';
          this.selection = { atoms: new Set(), bonds: new Set(), objects: new Set() };
          this.selection[key].add(d.hit.id);
          if (d.hit.type === 'bond') { const b = mol.getBond(d.hit.id); this.selection.atoms.add(b.a1); this.selection.atoms.add(b.a2); }
          this.render();
          this.dispatchEvent(new CustomEvent('selection'));
        }
        break;
      }
      case 'rotate': {
        if (d.moved) { this.undoStack.push(d.snap); this.redoStack = []; this._changed(); }
        break;
      }
      case 'rubber': {
        this.render();
        this.dispatchEvent(new CustomEvent('selection'));
        break;
      }
      case 'erase': {
        if (d.erased) { this.undoStack.push(d.snap); this.redoStack = []; this._changed(); }
        break;
      }
      case 'bond': this._finishBond(d, p); break;
      case 'atom': this._finishAtom(d, p); break;
      case 'chain': {
        if (!d.moved || !d.pts) { this._finishBond({ ...d, moved: false }, p, { order: 1, stereo: 'none' }); break; }
        const pts = d.pts;
        this.change(() => {
          let prev = d.hit?.type === 'atom' ? d.hit.id : this._atomAt(pts[0]).id;
          for (let i = 1; i < pts.length; i++) {
            const a = this._atomAt(pts[i]);
            if (a.id !== prev && !mol.bondBetween(prev, a.id)) mol.addBond(prev, a.id, 1);
            prev = a.id;
          }
        });
        break;
      }
      case 'ring': {
        const tmpl = this.tool.template;
        if (!tmpl) break;
        this.change(() => this.placeRing(d.hit, d.start, tmpl));
        this._updateHover(p);
        break;
      }
      case 'charge': {
        if (d.hit?.type === 'atom') this.change(() => this.changeCharge(d.hit.id, this.tool.delta || 1));
        break;
      }
      case 'text': {
        if (d.hit?.type === 'atom') this.editAtomLabel(d.hit.id);
        else if (d.hit?.type === 'object' && this.doc.getObject(d.hit.id).type === 'text') this.editTextObject(d.hit.id);
        else this.editTextObject(null, d.start);
        break;
      }
      case 'arrow': {
        const style = this.tool.style || 'forward';
        const end = d.moved && d.end ? d.end : add(d.start, { x: style === 'curved' || style === 'fishhook' ? 1.5 : 2.5, y: 0 });
        if (d.hit?.type === 'object' && !d.moved) {
          // Clicking a curved arrow flips its curvature.
          const o = this.doc.getObject(d.hit.id);
          if (o.type === 'arrow' && (o.style === 'curved' || o.style === 'fishhook')) { this.change(() => { o.curve = -(o.curve ?? 0.4); }); break; }
        }
        this.change(() => this.doc.addObject({ type: 'arrow', style, x1: d.start.x, y1: d.start.y, x2: end.x, y2: end.y, curve: 0.4 }));
        break;
      }
      default: break;
    }
    this.render();
  }

  _finishBond(d, p, toolOverride = null) {
    const mol = this.doc.mol;
    const tool = toolOverride ?? this.tool;
    const order = tool.order || 1, stereo = tool.stereo || 'none';
    this.change(() => {
      if (!d.moved) {
        if (d.hit?.type === 'bond') { if (!toolOverride) this._applyBondTool(mol.getBond(d.hit.id)); return; }
        if (d.hit?.type === 'atom') { this.addBondFromAtom(d.hit.id, order, stereo); return; }
        if (d.hit?.type === 'object') return;
        const a = mol.addAtom({ x: d.start.x, y: d.start.y });
        this.addBondFromAtom(a.id, order, stereo);
        return;
      }
      const fromId = d.hit?.type === 'atom' ? d.hit.id : this._atomAt(d.start).id;
      const end = d.end ?? this._bondDragEnd(mol.getAtom(fromId), p, false, fromId);
      const toId = end.atomId ?? mol.addAtom({ x: end.pos.x, y: end.pos.y }).id;
      if (toId === fromId) return;
      const existing = mol.bondBetween(fromId, toId);
      if (existing) { this._applyBondTool(existing); return; }
      mol.addBond(fromId, toId, order, stereo);
    });
  }

  _finishAtom(d, p) {
    const mol = this.doc.mol;
    const el = this.element;
    this.change(() => {
      if (!d.moved) {
        if (d.hit?.type === 'atom') {
          const a = mol.getAtom(d.hit.id);
          if (a.el === el && el !== 'C') this.addBondFromAtom(a.id, 1, 'none', el);
          else this.setAtomElement(a.id, el);
        } else if (!d.hit) {
          mol.addAtom({ el, x: d.start.x, y: d.start.y });
        }
        return;
      }
      const fromId = d.hit?.type === 'atom' ? d.hit.id : mol.addAtom({ x: d.start.x, y: d.start.y }).id;
      const end = d.end ?? this._bondDragEnd(mol.getAtom(fromId), p, false, fromId);
      const toId = end.atomId ?? mol.addAtom({ el, x: end.pos.x, y: end.pos.y }).id;
      if (toId !== fromId && !mol.bondBetween(fromId, toId)) mol.addBond(fromId, toId, 1);
    });
  }

  _onDblClick(e) {
    const p = this.toModel(e.clientX, e.clientY);
    const hit = this.hitTest(p);
    if (this.tool.name === 'select' && hit && hit.type !== 'object') {
      const atomId = hit.type === 'atom' ? hit.id : this.doc.mol.getBond(hit.id).a1;
      const ids = this.doc.mol.connectedAtoms(atomId);
      this.selection.atoms = new Set(ids);
      this._selectBondsBetweenSelectedAtoms();
      this.render();
      this.dispatchEvent(new CustomEvent('selection'));
    } else if (this.tool.name === 'select' && hit?.type === 'object' && this.doc.getObject(hit.id).type === 'text') {
      this.editTextObject(hit.id);
    }
  }

  // ========================================================= text editing
  _openTextEditor(screen, value, onCommit, { multiline = false } = {}) {
    this._closeTextEditor(false);
    const host = this.svg.parentElement;
    const el = document.createElement(multiline ? 'textarea' : 'input');
    el.className = 'bz-inline-editor';
    el.value = value;
    el.style.left = `${screen.x}px`;
    el.style.top = `${screen.y}px`;
    host.appendChild(el);
    el.focus();
    el.select();
    const commit = (ok) => {
      if (!this._textEditor) return;
      const v = el.value;
      this._textEditor = null;
      el.remove();
      if (ok) onCommit(v);
    };
    el.addEventListener('keydown', (ev) => {
      ev.stopPropagation();
      if (ev.key === 'Enter' && !(multiline && ev.shiftKey)) { ev.preventDefault(); commit(true); }
      else if (ev.key === 'Escape') { ev.preventDefault(); commit(false); }
    });
    el.addEventListener('blur', () => commit(true));
    this._textEditor = { el, commit };
  }

  _closeTextEditor(ok = true) {
    if (this._textEditor) this._textEditor.commit(ok);
  }

  editAtomLabel(atomId) {
    const a = this.doc.mol.getAtom(atomId);
    const s = this.toScreen(a);
    const r = this.svg.getBoundingClientRect();
    const host = this.svg.parentElement.getBoundingClientRect();
    this._openTextEditor({ x: s.x + r.left - host.left - 14, y: s.y + r.top - host.top - 12 }, a.label || (a.el === 'C' ? '' : a.el), (v) => {
      if (v.trim()) this.change(() => this.applyAtomLabel(atomId, v));
    });
  }

  editTextObject(objId, p = null) {
    const o = objId ? this.doc.getObject(objId) : null;
    const pos = o ? { x: o.x, y: o.y } : p;
    const s = this.toScreen(pos);
    const r = this.svg.getBoundingClientRect();
    const host = this.svg.parentElement.getBoundingClientRect();
    this._openTextEditor({ x: s.x + r.left - host.left, y: s.y + r.top - host.top - 18 }, o ? o.text : '', (v) => {
      this.change(() => {
        if (o) { if (v.trim()) o.text = v; else this.doc.removeObject(o.id); }
        else if (v.trim()) this.doc.addObject({ type: 'text', x: pos.x, y: pos.y, text: v, size: 1 });
      });
    }, { multiline: true });
  }

  // ============================================================== hotkeys
  /** Handle a keydown event. Returns true if handled. */
  handleKey(e) {
    if (this._textEditor) return false;
    const k = e.key;
    const mod = e.ctrlKey || e.metaKey;
    const mol = this.doc.mol;
    if (k === ' ') { this.spaceDown = true; this.svg.style.cursor = 'grab'; return true; }
    if (mod) return false; // handled by the app (menus)
    const h = this.hover;
    if (k === 'Delete' || k === 'Backspace') {
      if (this.hasSelection()) this.deleteSelection();
      else if (h) this.change(() => { this._deleteHit(h); this.hover = null; });
      return true;
    }
    if (k === 'Escape') { this.clearSelection(); this.setTool({ name: 'select' }); return true; }
    if (k.startsWith('Arrow') && this.hasSelection()) {
      const step = e.shiftKey ? 0.5 : 0.1;
      const dx = k === 'ArrowLeft' ? -step : k === 'ArrowRight' ? step : 0;
      const dy = k === 'ArrowUp' ? -step : k === 'ArrowDown' ? step : 0;
      this.nudge(dx, dy);
      return true;
    }
    if (h?.type === 'atom') {
      const id = h.id;
      if (k === '+' || k === '=') { this.change(() => this.changeCharge(id, 1)); return true; }
      if (k === '-' || k === '_') { this.change(() => this.changeCharge(id, -1)); return true; }
      if (k === '1' || k === '2' || k === '3') { this.change(() => this.addBondFromAtom(id, Number(k))); return true; }
      if (k === 'Enter') { this.editAtomLabel(id); return true; }
      const el = ELEMENT_KEYS[k];
      if (el) { this.change(() => this.setAtomElement(id, el)); return true; }
      return false;
    }
    if (h?.type === 'bond') {
      const b = mol.getBond(h.id);
      if (k === '1' || k === '2' || k === '3') { this.change(() => { b.order = Number(k); b.stereo = 'none'; }); return true; }
      if (k === 'w') { this.change(() => { if (b.stereo === 'wedge') [b.a1, b.a2] = [b.a2, b.a1]; b.order = 1; b.stereo = 'wedge'; mol.invalidate(); }); return true; }
      if (k === 'h' || k === 'q') { this.change(() => { if (b.stereo === 'hash') [b.a1, b.a2] = [b.a2, b.a1]; b.order = 1; b.stereo = 'hash'; mol.invalidate(); }); return true; }
      if (k === 'y') { this.change(() => { b.order = 1; b.stereo = 'wavy'; }); return true; }
      if (k === '0') { this.change(() => { b.stereo = 'none'; }); return true; }
      return false;
    }
    return false;
  }

  handleKeyUp(e) {
    if (e.key === ' ') { this.spaceDown = false; this.svg.style.cursor = this._cursorFor(this.hover); }
  }
}

export { freeDirection };
