// Core molecular graph model.
//
// Conventions:
//  * Coordinates are in "bond length units": a standard bond is 1.0 long.
//  * The y axis points DOWN (screen convention). File formats with y-up
//    (MOL files) must flip y on import/export.
//  * Bonds are stored in Kekulé form (order 1, 2 or 3). Aromaticity is a
//    derived property (see rings.js).
//  * atom.hCount === null means "compute implicit hydrogens from valence";
//    a number forces that many hydrogens (e.g. from SMILES brackets).

import { VALENCES, ELEMENTS } from './elements.js';

export const BOND_STEREO = Object.freeze({
  NONE: 'none', WEDGE: 'wedge', HASH: 'hash', WAVY: 'wavy', BOLD: 'bold', DASHED: 'dashed',
});

export class Molecule {
  constructor() {
    /** @type {Map<number, Atom>} */
    this.atoms = new Map();
    /** @type {Map<number, Bond>} */
    this.bonds = new Map();
    this.nextId = 1;
    this._adj = null;
  }

  // ---------------------------------------------------------------- mutation
  addAtom(props = {}) {
    const atom = {
      id: props.id ?? this.nextId++,
      el: props.el ?? 'C',
      x: props.x ?? 0,
      y: props.y ?? 0,
      charge: props.charge ?? 0,
      isotope: props.isotope ?? null,
      hCount: props.hCount ?? null,
      radical: props.radical ?? 0, // number of unpaired electrons (0, 1, 2)
      label: props.label ?? null, // optional display label / abbreviation (e.g. "OMe")
      mapNo: props.mapNo ?? 0,
    };
    if (atom.id >= this.nextId) this.nextId = atom.id + 1;
    this.atoms.set(atom.id, atom);
    this._adj = null;
    return atom;
  }

  addBond(a1, a2, order = 1, stereo = BOND_STEREO.NONE, props = {}) {
    if (a1 === a2) throw new Error('Cannot bond an atom to itself');
    if (!this.atoms.has(a1) || !this.atoms.has(a2)) throw new Error('Bond references unknown atom');
    const existing = this.bondBetween(a1, a2);
    if (existing) return existing;
    const bond = {
      id: props.id ?? this.nextId++,
      a1, a2, order, stereo,
      // Double bond placement: 'auto' | 'center' | 'left' | 'right'
      position: props.position ?? 'auto',
    };
    if (bond.id >= this.nextId) this.nextId = bond.id + 1;
    this.bonds.set(bond.id, bond);
    this._adj = null;
    return bond;
  }

  removeAtom(id) {
    for (const b of this.bondsOf(id)) this.bonds.delete(b.id);
    this.atoms.delete(id);
    this._adj = null;
  }

  removeBond(id) {
    this.bonds.delete(id);
    this._adj = null;
  }

  /** Call after mutating bond endpoints directly. */
  invalidate() {
    this._adj = null;
  }

  // ----------------------------------------------------------------- queries
  getAtom(id) { return this.atoms.get(id); }
  getBond(id) { return this.bonds.get(id); }
  atomList() { return [...this.atoms.values()]; }
  bondList() { return [...this.bonds.values()]; }
  get atomCount() { return this.atoms.size; }
  get bondCount() { return this.bonds.size; }
  isEmpty() { return this.atoms.size === 0; }

  _buildAdj() {
    const adj = new Map();
    for (const id of this.atoms.keys()) adj.set(id, []);
    for (const b of this.bonds.values()) {
      adj.get(b.a1)?.push(b);
      adj.get(b.a2)?.push(b);
    }
    this._adj = adj;
  }

  /** Bonds incident on an atom. */
  bondsOf(atomId) {
    if (!this._adj) this._buildAdj();
    return this._adj.get(atomId) || [];
  }

  /** Neighbor atom ids. */
  neighbors(atomId) {
    return this.bondsOf(atomId).map((b) => (b.a1 === atomId ? b.a2 : b.a1));
  }

  degree(atomId) { return this.bondsOf(atomId).length; }

  bondBetween(a, b) {
    for (const bond of this.bondsOf(a)) {
      if ((bond.a1 === a && bond.a2 === b) || (bond.a1 === b && bond.a2 === a)) return bond;
    }
    return null;
  }

  otherAtom(bond, atomId) { return bond.a1 === atomId ? bond.a2 : bond.a1; }

  /** Sum of bond orders to explicit neighbours. */
  bondOrderSum(atomId) {
    let s = 0;
    for (const b of this.bondsOf(atomId)) s += b.order;
    return s;
  }

  /** Number of implicit (non-drawn) hydrogens on an atom. */
  implicitH(atomId) {
    const atom = this.atoms.get(atomId);
    if (!atom) return 0;
    if (atom.hCount !== null && atom.hCount !== undefined) return atom.hCount;
    return computeImplicitH(atom.el, atom.charge, this.bondOrderSum(atomId), atom.radical);
  }

  /** Total H count = implicit H + explicitly drawn H neighbours. */
  totalH(atomId) {
    let n = this.implicitH(atomId);
    for (const nb of this.neighbors(atomId)) if (this.atoms.get(nb).el === 'H') n++;
    return n;
  }

  /** Returns true when an atom exceeds every allowed valence. */
  hasValenceError(atomId) {
    const atom = this.atoms.get(atomId);
    const vals = VALENCES[atom.el];
    if (!vals) return false;
    const used = this.bondOrderSum(atomId) + (atom.hCount ?? 0) + atom.radical;
    const max = Math.max(...vals.map((v) => adjustValence(atom.el, v, atom.charge)));
    return used > max;
  }

  /** Connected components as arrays of atom ids. */
  fragments() {
    const seen = new Set();
    const out = [];
    for (const id of this.atoms.keys()) {
      if (seen.has(id)) continue;
      const comp = [];
      const stack = [id];
      seen.add(id);
      while (stack.length) {
        const cur = stack.pop();
        comp.push(cur);
        for (const nb of this.neighbors(cur)) {
          if (!seen.has(nb)) { seen.add(nb); stack.push(nb); }
        }
      }
      out.push(comp);
    }
    return out;
  }

  /** Atom ids connected to `startId` (including it). */
  connectedAtoms(startId) {
    const seen = new Set([startId]);
    const stack = [startId];
    while (stack.length) {
      const cur = stack.pop();
      for (const nb of this.neighbors(cur)) if (!seen.has(nb)) { seen.add(nb); stack.push(nb); }
    }
    return seen;
  }

  bbox(atomIds = null) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    const ids = atomIds ?? this.atoms.keys();
    for (const id of ids) {
      const a = this.atoms.get(id);
      if (!a) continue;
      minX = Math.min(minX, a.x); maxX = Math.max(maxX, a.x);
      minY = Math.min(minY, a.y); maxY = Math.max(maxY, a.y);
    }
    if (minX === Infinity) return null;
    return { minX, minY, maxX, maxY, width: maxX - minX, height: maxY - minY, cx: (minX + maxX) / 2, cy: (minY + maxY) / 2 };
  }

  /** Average bond length (1.0 when there are no bonds). */
  averageBondLength() {
    if (this.bonds.size === 0) return 1;
    let s = 0;
    for (const b of this.bonds.values()) {
      const p = this.atoms.get(b.a1), q = this.atoms.get(b.a2);
      s += Math.hypot(p.x - q.x, p.y - q.y);
    }
    return s / this.bonds.size || 1;
  }

  translate(dx, dy, atomIds = null) {
    for (const id of atomIds ?? this.atoms.keys()) {
      const a = this.atoms.get(id);
      a.x += dx; a.y += dy;
    }
  }

  scale(f, atomIds = null) {
    for (const id of atomIds ?? this.atoms.keys()) {
      const a = this.atoms.get(id);
      a.x *= f; a.y *= f;
    }
  }

  // ------------------------------------------------------------ composition
  clone() {
    const m = new Molecule();
    for (const a of this.atoms.values()) m.atoms.set(a.id, { ...a });
    for (const b of this.bonds.values()) m.bonds.set(b.id, { ...b });
    m.nextId = this.nextId;
    return m;
  }

  /** Extract a sub-molecule containing only the given atom ids (ids preserved). */
  subMolecule(atomIds) {
    const set = new Set(atomIds);
    const m = new Molecule();
    for (const id of set) { const a = this.atoms.get(id); if (a) m.atoms.set(id, { ...a }); }
    for (const b of this.bonds.values()) if (set.has(b.a1) && set.has(b.a2)) m.bonds.set(b.id, { ...b });
    m.nextId = this.nextId;
    return m;
  }

  /** Merge another molecule into this one; returns Map(oldId -> newId) for atoms. */
  merge(other) {
    const map = new Map();
    for (const a of other.atoms.values()) {
      const { id, ...rest } = a;
      map.set(id, this.addAtom(rest).id);
    }
    for (const b of other.bonds.values()) {
      this.addBond(map.get(b.a1), map.get(b.a2), b.order, b.stereo, { position: b.position });
    }
    return map;
  }

  toJSON() {
    return {
      atoms: this.atomList().map((a) => ({ ...a })),
      bonds: this.bondList().map((b) => ({ ...b })),
    };
  }

  static fromJSON(json) {
    const m = new Molecule();
    for (const a of json.atoms || []) m.addAtom(a);
    for (const b of json.bonds || []) m.addBond(b.a1, b.a2, b.order, b.stereo, b);
    return m;
  }
}

/** Valence adjusted for formal charge. */
export function adjustValence(el, v, charge) {
  const info = ELEMENTS[el];
  if (!info) return v;
  const col = info.col;
  if (col === 13) return v - charge; // B- -> 4
  if (col === 14) return v - Math.abs(charge); // C+ / C- -> 3
  if (col >= 15 && col <= 17) return v + charge; // N+ -> 4, O- -> 1
  return v;
}

export function computeImplicitH(el, charge, bondOrderSum, radical = 0) {
  const vals = VALENCES[el];
  if (!vals) return 0;
  const used = bondOrderSum + radical;
  for (const v of vals) {
    const av = adjustValence(el, v, charge);
    if (av >= used) return av - used;
  }
  return 0;
}
