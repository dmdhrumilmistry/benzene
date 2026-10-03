// MDL Molfile (V2000 read/write, basic V3000 read) and SDF (first record).

import { Molecule, BOND_STEREO } from './molecule.js';
import { isElement, ELEMENTS } from './elements.js';
import { kekulize } from './smiles.js';
import { layoutMolecule } from './layout.js';

const EXPORT_BOND_LENGTH = 1.5;
const CHARGE_CODES = { 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3 };

/**
 * Parse a V2000/V3000 molfile (or the first record of an SDF) into a Molecule.
 * Coordinates are converted to y-down with median bond length 1.0; molecules
 * without 2D coordinates get an automatic layout.
 * @param {string} text
 * @returns {Molecule}
 * @throws {Error} on malformed input
 */
export function readMolfile(text) {
  let src = String(text ?? '').replace(/\r\n?/g, '\n');
  const sdfEnd = src.indexOf('\n$$$$');
  if (sdfEnd >= 0) src = src.slice(0, sdfEnd);
  const lines = src.split('\n');
  if (lines.length < 4) throw new Error('Invalid molfile: too few lines');
  const counts = lines[3];
  const mol = new Molecule();
  const aromaticBonds = [];
  const valenceFields = new Map();
  if (/V3000/i.test(counts)) readV3000(lines, mol, aromaticBonds);
  else readV2000(lines, mol, aromaticBonds, valenceFields);

  if (aromaticBonds.length) kekulize(mol, aromaticBonds);
  // Explicit valence field: fixes the hydrogen count.
  for (const [id, v] of valenceFields) {
    const atom = mol.getAtom(id);
    atom.hCount = v === 15 ? 0 : Math.max(0, v - mol.bondOrderSum(id));
  }
  normalizeCoordinates(mol);
  return mol;
}

function int(s, def = 0) {
  const n = parseInt(s, 10);
  return Number.isFinite(n) ? n : def;
}

function makeAtom(mol, sym, x, y, extra = {}) {
  const props = { x, y: -y, ...extra };
  if (sym === 'D' || sym === 'T') {
    props.el = 'H';
    props.isotope = sym === 'D' ? 2 : 3;
  } else if (isElement(sym)) {
    props.el = sym;
  } else {
    // Query / pseudo atoms (R#, A, Q, *, L...) are kept as labelled carbons.
    props.el = 'C';
    props.label = sym;
    props.hCount = 0;
  }
  return mol.addAtom(props);
}

function applyStereo(code, isV3000) {
  if (isV3000) {
    if (code === 1) return BOND_STEREO.WEDGE;
    if (code === 3) return BOND_STEREO.HASH;
    if (code === 2) return BOND_STEREO.WAVY;
    return BOND_STEREO.NONE;
  }
  if (code === 1) return BOND_STEREO.WEDGE;
  if (code === 6) return BOND_STEREO.HASH;
  if (code === 4) return BOND_STEREO.WAVY;
  return BOND_STEREO.NONE;
}

function addMolBond(mol, a1, a2, type, stereo, aromaticBonds) {
  if (!a1 || !a2 || a1 === a2) throw new Error('Invalid molfile: bond references unknown atom');
  const order = type >= 1 && type <= 3 ? type : 1;
  const bond = mol.addBond(a1, a2, order, stereo);
  if (type === 4) aromaticBonds.push(bond.id);
  return bond;
}

function readV2000(lines, mol, aromaticBonds, valenceFields) {
  const counts = lines[3];
  let nAtoms = int(counts.slice(0, 3), NaN);
  let nBonds = int(counts.slice(3, 6), NaN);
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds)) {
    const t = counts.trim().split(/\s+/);
    nAtoms = int(t[0], NaN); nBonds = int(t[1], NaN);
  }
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds)) throw new Error('Invalid molfile: bad counts line');
  if (lines.length < 4 + nAtoms + nBonds) throw new Error('Invalid molfile: truncated atom/bond block');
  const ids = [];
  for (let k = 0; k < nAtoms; k++) {
    const line = lines[4 + k];
    let x, y, sym, dd, ccc, vvv;
    if (line.length >= 34) {
      x = parseFloat(line.slice(0, 10)); y = parseFloat(line.slice(10, 20));
      sym = line.slice(31, 34).trim();
      dd = int(line.slice(34, 36)); ccc = int(line.slice(36, 39)); vvv = int(line.slice(48, 51));
    }
    if (!sym || !Number.isFinite(x) || !Number.isFinite(y)) {
      const t = line.trim().split(/\s+/);
      x = parseFloat(t[0]); y = parseFloat(t[1]); sym = t[3];
      dd = int(t[4]); ccc = int(t[5]); vvv = int(t[9]);
    }
    if (!sym || !Number.isFinite(x) || !Number.isFinite(y)) throw new Error(`Invalid molfile: bad atom line ${5 + k}`);
    const extra = {};
    if (CHARGE_CODES[ccc]) extra.charge = CHARGE_CODES[ccc];
    if (ccc === 4) extra.radical = 1;
    const atom = makeAtom(mol, sym, x, y, extra);
    if (dd && ELEMENTS[atom.el] && !atom.label) atom.isotope = Math.round(ELEMENTS[atom.el].mass) + dd;
    if (vvv) valenceFields.set(atom.id, vvv);
    ids.push(atom.id);
  }
  for (let k = 0; k < nBonds; k++) {
    const line = lines[4 + nAtoms + k];
    let a1 = int(line.slice(0, 3), NaN), a2 = int(line.slice(3, 6), NaN);
    let type = int(line.slice(6, 9), NaN), st = int(line.slice(9, 12));
    if (!Number.isFinite(a1) || !Number.isFinite(a2) || !Number.isFinite(type)) {
      const t = line.trim().split(/\s+/);
      a1 = int(t[0]); a2 = int(t[1]); type = int(t[2]); st = int(t[3]);
    }
    addMolBond(mol, ids[a1 - 1], ids[a2 - 1], type, applyStereo(st, false), aromaticBonds);
  }
  // Properties block.
  let resetDone = false;
  const resetChargeRad = () => {
    if (resetDone) return;
    resetDone = true;
    for (const a of mol.atoms.values()) { a.charge = 0; a.radical = 0; }
  };
  for (let k = 4 + nAtoms + nBonds; k < lines.length; k++) {
    const line = lines[k];
    if (line.startsWith('M  END')) break;
    if (line.startsWith('A  ')) { k++; continue; }
    const tag = line.slice(0, 6);
    if (tag !== 'M  CHG' && tag !== 'M  ISO' && tag !== 'M  RAD') continue;
    const t = line.slice(6).trim().split(/\s+/).map((v) => int(v));
    const n = t[0];
    if (tag !== 'M  ISO') resetChargeRad();
    for (let j = 0; j < n; j++) {
      const atom = mol.getAtom(ids[t[1 + 2 * j] - 1]);
      const v = t[2 + 2 * j];
      if (!atom) continue;
      if (tag === 'M  CHG') atom.charge = v;
      else if (tag === 'M  ISO') atom.isotope = v;
      else atom.radical = v === 2 ? 1 : v === 1 || v === 3 ? 2 : 0;
    }
  }
}

function readV3000(lines, mol, aromaticBonds) {
  // Join continuation lines (ending in '-').
  const body = [];
  for (let k = 4; k < lines.length; k++) {
    let line = lines[k];
    if (line.startsWith('M  END')) break;
    if (!line.startsWith('M  V30 ')) continue;
    line = line.slice(7);
    while (line.endsWith('-') && k + 1 < lines.length) line = line.slice(0, -1) + lines[++k].replace(/^M {2}V30 /, '');
    body.push(line.trim());
  }
  const idx = new Map();
  let section = null;
  for (const line of body) {
    if (/^BEGIN\s+(\w+)/.test(line)) { section = line.split(/\s+/)[1]; continue; }
    if (/^END\s+/.test(line)) { section = null; continue; }
    const tokens = line.match(/(?:[^\s"(]+(?:\([^)]*\))?|"[^"]*")+/g) || [];
    const kv = {};
    for (const tok of tokens.slice(1)) {
      const m = /^([A-Z]+)=(.*)$/.exec(tok);
      if (m) kv[m[1]] = m[2];
    }
    if (section === 'ATOM') {
      const [n, sym, x, y] = tokens;
      const extra = {};
      if (kv.CHG) extra.charge = int(kv.CHG);
      if (kv.MASS) extra.isotope = int(kv.MASS);
      if (kv.RAD) { const r = int(kv.RAD); extra.radical = r === 2 ? 1 : r ? 2 : 0; }
      const atom = makeAtom(mol, sym, parseFloat(x), parseFloat(y), extra);
      if (kv.VAL) {
        const v = int(kv.VAL);
        atom.hCount = v === -1 ? 0 : null;
        if (v > 0) atom._val = v;
      }
      idx.set(int(n), atom.id);
    } else if (section === 'BOND') {
      const [, type, a1, a2] = tokens;
      addMolBond(mol, idx.get(int(a1)), idx.get(int(a2)), int(type), applyStereo(int(kv.CFG), true), aromaticBonds);
    }
  }
  for (const atom of mol.atoms.values()) {
    if (atom._val) { atom.hCount = Math.max(0, atom._val - mol.bondOrderSum(atom.id)); delete atom._val; }
  }
}

function normalizeCoordinates(mol) {
  const atoms = mol.atomList();
  if (!atoms.length) return;
  const allSame = atoms.every((a) => Math.abs(a.x - atoms[0].x) < 1e-6 && Math.abs(a.y - atoms[0].y) < 1e-6);
  if (allSame && atoms.length > 1) {
    layoutMolecule(mol);
    return;
  }
  const lens = mol.bondList()
    .map((b) => { const p = mol.getAtom(b.a1), q = mol.getAtom(b.a2); return Math.hypot(p.x - q.x, p.y - q.y); })
    .filter((l) => l > 1e-6)
    .sort((a, b) => a - b);
  if (lens.length) {
    const med = lens[Math.floor(lens.length / 2)];
    mol.scale(1 / med);
  }
  const bb = mol.bbox();
  mol.translate(-bb.cx, -bb.cy);
}

// ------------------------------------------------------------------- writer

const pad = (v, w) => String(v).padStart(w);
const fmt = (v) => (Math.abs(v) < 5e-5 ? 0 : v).toFixed(4).padStart(10);

/**
 * Write a V2000 molfile. Coordinates are y-flipped and scaled to 1.5 Å bonds.
 * Charges, isotopes and radicals are written to M  CHG / M  ISO / M  RAD.
 * @param {Molecule} mol
 * @param {{name?: string}} [opts]
 * @returns {string} molfile text ending with 'M  END'
 */
export function writeMolfile(mol, { name = '' } = {}) {
  const atoms = mol.atomList();
  const bonds = mol.bondList();
  if (atoms.length > 999 || bonds.length > 999) throw new Error('V2000 molfiles support at most 999 atoms and bonds');
  const index = new Map(atoms.map((a, k) => [a.id, k + 1]));
  const out = [];
  out.push(String(name).replace(/[\r\n]/g, ' ').slice(0, 80));
  out.push('  Benzene           2D');
  out.push('');
  out.push(`${pad(atoms.length, 3)}${pad(bonds.length, 3)}  0  0  0  0  0  0  0  0999 V2000`);
  const chg = [], iso = [], rad = [];
  for (const a of atoms) {
    const i = index.get(a.id);
    const code = { 3: 1, 2: 2, 1: 3, [-1]: 5, [-2]: 6, [-3]: 7 }[a.charge] ?? 0;
    const sym = (a.label && !isElement(a.label) && a.hCount === 0 && a.el === 'C') ? 'R' : a.el;
    let vvv = 0;
    if (a.hCount !== null && a.hCount !== undefined) {
      const auto = computeAutoH(mol, a);
      if (auto !== a.hCount) {
        const v = mol.bondOrderSum(a.id) + a.hCount;
        vvv = v === 0 ? 15 : Math.min(14, v);
      }
    }
    out.push(`${fmt(a.x * EXPORT_BOND_LENGTH)}${fmt(-a.y * EXPORT_BOND_LENGTH)}${fmt(0)} ${sym.padEnd(3)} 0${pad(code, 3)}  0  0  0${pad(vvv, 3)}  0  0  0  0  0  0`);
    if (a.charge) chg.push([i, a.charge]);
    if (a.isotope) iso.push([i, a.isotope]);
    if (a.radical) rad.push([i, a.radical === 1 ? 2 : 3]);
  }
  for (const b of bonds) {
    const st = b.stereo === BOND_STEREO.WEDGE ? 1 : b.stereo === BOND_STEREO.HASH ? 6 : b.stereo === BOND_STEREO.WAVY ? 4 : 0;
    const type = b.order >= 1 && b.order <= 3 ? b.order : 1;
    out.push(`${pad(index.get(b.a1), 3)}${pad(index.get(b.a2), 3)}${pad(type, 3)}${pad(st, 3)}  0  0  0`);
  }
  const prop = (tag, list) => {
    for (let k = 0; k < list.length; k += 8) {
      const chunk = list.slice(k, k + 8);
      out.push(`M  ${tag}${pad(chunk.length, 3)}${chunk.map(([i, v]) => ` ${pad(i, 3)} ${pad(v, 3)}`).join('')}`);
    }
  };
  prop('CHG', chg);
  prop('ISO', iso);
  prop('RAD', rad);
  out.push('M  END');
  return out.join('\n');
}

function computeAutoH(mol, atom) {
  const saved = atom.hCount;
  atom.hCount = null;
  const h = mol.implicitH(atom.id);
  atom.hCount = saved;
  return h;
}
