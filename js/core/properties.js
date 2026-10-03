// Molecular property calculations: formula, masses, elemental analysis,
// isotope pattern, unsaturation, H-bond donors/acceptors, rotatable bonds, TPSA.

import { ELEMENTS } from './elements.js';
import { analyzeRings } from './rings.js';

const ELECTRON_MASS = 0.00054858;

/** Element counts (including implicit hydrogens) for the given atoms. */
export function elementCounts(mol, atomIds = null) {
  const counts = {};
  for (const id of atomIds ?? mol.atoms.keys()) {
    const a = mol.getAtom(id);
    if (!ELEMENTS[a.el]) continue;
    counts[a.el] = (counts[a.el] || 0) + 1;
    const h = mol.implicitH(id);
    if (h) counts.H = (counts.H || 0) + h;
  }
  return counts;
}

export function totalCharge(mol, atomIds = null) {
  let c = 0;
  for (const id of atomIds ?? mol.atoms.keys()) c += mol.getAtom(id).charge || 0;
  return c;
}

/** Hill-order list of [symbol, count]. */
export function hillOrder(counts) {
  const syms = Object.keys(counts);
  const hasC = 'C' in counts;
  syms.sort((a, b) => {
    if (hasC) {
      const rank = (s) => (s === 'C' ? 0 : s === 'H' ? 1 : 2);
      if (rank(a) !== rank(b)) return rank(a) - rank(b);
    }
    return a < b ? -1 : a > b ? 1 : 0;
  });
  return syms.map((s) => [s, counts[s]]);
}

/** Plain-text formula, e.g. "C9H8O4" or "C2H3O2-". */
export function formulaString(counts, charge = 0) {
  let s = hillOrder(counts).map(([el, n]) => el + (n > 1 ? n : '')).join('');
  if (charge) s += (Math.abs(charge) > 1 ? Math.abs(charge) : '') + (charge > 0 ? '+' : '-');
  return s;
}

/** HTML formula with subscripts / superscripted charge. */
export function formulaHTML(counts, charge = 0) {
  let s = hillOrder(counts).map(([el, n]) => el + (n > 1 ? `<sub>${n}</sub>` : '')).join('');
  if (charge) s += `<sup>${Math.abs(charge) > 1 ? Math.abs(charge) : ''}${charge > 0 ? '+' : '−'}</sup>`;
  return s;
}

export function molecularWeight(counts) {
  let m = 0;
  for (const [el, n] of Object.entries(counts)) m += ELEMENTS[el].mass * n;
  return m;
}

export function exactMass(counts, charge = 0) {
  let m = 0;
  for (const [el, n] of Object.entries(counts)) m += ELEMENTS[el].monoisotopic * n;
  return m - charge * ELECTRON_MASS;
}

/** Mass percent of each element (Hill order). */
export function elementalAnalysis(counts) {
  const mw = molecularWeight(counts);
  if (!mw) return [];
  return hillOrder(counts).map(([el, n]) => ({ el, percent: (ELEMENTS[el].mass * n * 100) / mw }));
}

/**
 * Isotope distribution. Returns peaks [{mass, abundance}] with the most
 * intense peak normalised to 100, merged to `resolution` Da and pruned below
 * `threshold` (% of max).
 */
export function isotopePattern(counts, charge = 0, { threshold = 0.05, resolution = 0.01, maxPeaks = 40 } = {}) {
  let dist = [[0, 1]];
  const prune = (d) => {
    d.sort((a, b) => a[0] - b[0]);
    const merged = [];
    for (const [m, p] of d) {
      const last = merged[merged.length - 1];
      if (last && m - last[0] < resolution) {
        const tot = last[1] + p;
        last[0] = (last[0] * last[1] + m * p) / tot;
        last[1] = tot;
      } else merged.push([m, p]);
    }
    const max = Math.max(...merged.map((x) => x[1]));
    return merged.filter((x) => x[1] >= max * 1e-6);
  };
  for (const [el, n] of Object.entries(counts)) {
    const info = ELEMENTS[el];
    const iso = info.isotopes || [[info.monoisotopic, 1]];
    for (let i = 0; i < n; i++) {
      const next = [];
      for (const [m1, p1] of dist) for (const [m2, p2] of iso) next.push([m1 + m2, p1 * p2]);
      dist = prune(next);
    }
  }
  const z = Math.abs(charge) || 1;
  const max = Math.max(...dist.map((x) => x[1]));
  return dist
    .map(([m, p]) => ({ mass: (m - charge * ELECTRON_MASS) / z, abundance: (p / max) * 100 }))
    .filter((x) => x.abundance >= threshold)
    .sort((a, b) => b.abundance - a.abundance)
    .slice(0, maxPeaks)
    .sort((a, b) => a.mass - b.mass);
}

/** Degree of unsaturation (rings + pi bonds) from formula. */
export function degreeOfUnsaturation(counts) {
  let c = 0, h = 0, n = 0;
  for (const [el, k] of Object.entries(counts)) {
    const v = ELEMENTS[el].valences?.[0];
    if (v === 4) c += k;
    else if (v === 1) h += k;
    else if (v === 3) n += k;
  }
  return c - h / 2 + n / 2 + 1;
}

function isHeteroWithH(mol, id) {
  const a = mol.getAtom(id);
  return (a.el === 'N' || a.el === 'O') && mol.totalH(id) > 0;
}

/** Lipinski-style H-bond donors: N-H and O-H groups (counted per atom). */
export function hBondDonors(mol, atomIds = null) {
  let n = 0;
  for (const id of atomIds ?? mol.atoms.keys()) if (isHeteroWithH(mol, id)) n++;
  return n;
}

/** Lipinski-style H-bond acceptors: N and O atoms. */
export function hBondAcceptors(mol, atomIds = null) {
  let n = 0;
  for (const id of atomIds ?? mol.atoms.keys()) {
    const el = mol.getAtom(id).el;
    if (el === 'N' || el === 'O') n++;
  }
  return n;
}

/** Rotatable bonds: non-ring single bonds between two non-terminal heavy atoms (excluding amide C-N and bonds to CX3/C≡). */
export function rotatableBonds(mol, ringInfo = analyzeRings(mol), atomIds = null) {
  const subset = atomIds ? new Set(atomIds) : null;
  const heavyDegree = (id) => mol.neighbors(id).filter((nb) => mol.getAtom(nb).el !== 'H').length;
  const isTripleAtom = (id) => mol.bondsOf(id).some((b) => b.order === 3);
  const isAmideCN = (b) => {
    const [c, n] = mol.getAtom(b.a1).el === 'C' ? [b.a1, b.a2] : [b.a2, b.a1];
    if (mol.getAtom(c).el !== 'C' || mol.getAtom(n).el !== 'N') return false;
    return mol.bondsOf(c).some((x) => x.order === 2 && mol.getAtom(mol.otherAtom(x, c)).el === 'O');
  };
  let count = 0;
  for (const b of mol.bonds.values()) {
    if (subset && !subset.has(b.a1)) continue;
    if (b.order !== 1 || ringInfo.bondRings.has(b.id)) continue;
    if (mol.getAtom(b.a1).el === 'H' || mol.getAtom(b.a2).el === 'H') continue;
    if (heavyDegree(b.a1) < 2 || heavyDegree(b.a2) < 2) continue;
    if (isTripleAtom(b.a1) || isTripleAtom(b.a2)) continue;
    if (isAmideCN(b)) continue;
    count++;
  }
  return count;
}

/**
 * Topological polar surface area (Ertl et al. 2000), N and O contributions only
 * (the commonly reported "TPSA" convention).
 */
export function tpsa(mol, ringInfo = analyzeRings(mol), atomIds = null) {
  const arom = ringInfo.aromatic;
  let total = 0;
  for (const id of atomIds ?? mol.atoms.keys()) {
    const a = mol.getAtom(id);
    if (a.el !== 'N' && a.el !== 'O') continue;
    const bonds = mol.bondsOf(id).filter((b) => mol.getAtom(mol.otherAtom(b, id)).el !== 'H');
    const h = mol.totalH(id);
    const isArom = arom.atoms.has(id);
    const nSingle = bonds.filter((b) => b.order === 1 && !arom.bonds.has(b.id)).length;
    const nDouble = bonds.filter((b) => b.order === 2 && !arom.bonds.has(b.id)).length;
    const nTriple = bonds.filter((b) => b.order === 3).length;
    const nArom = bonds.filter((b) => arom.bonds.has(b.id)).length;
    const in3Ring = (ringInfo.atomRings.get(id) || []).some((r) => ringInfo.rings[r].length === 3);
    const q = a.charge;
    let v = 0;
    if (a.el === 'N') {
      if (isArom) {
        if (q === 0 && h === 0 && nArom === 2 && nSingle === 0) v = 12.89;
        else if (q === 0 && h === 0 && nArom === 3) v = 4.41;
        else if (q === 0 && h === 0 && nArom === 2 && nSingle === 1) v = 4.93;
        else if (q === 0 && h === 0 && nArom === 2 && nDouble === 1) v = 8.39;
        else if (q === 0 && h === 1) v = 15.79;
        else if (q === 1 && h === 0 && nArom === 3) v = 4.10;
        else if (q === 1 && h === 0 && nArom === 2 && nSingle === 1) v = 3.88;
        else if (q === 1 && h === 1) v = 14.14;
        else v = 4.41;
      } else if (q === 0) {
        if (h === 0 && nSingle === 3) v = in3Ring ? 3.01 : 3.24;
        else if (h === 0 && nSingle === 1 && nDouble === 1) v = 12.36;
        else if (h === 0 && nTriple === 1) v = 23.79;
        else if (h === 0 && nSingle === 1 && nDouble === 2) v = 11.68;
        else if (h === 0 && nDouble === 1 && nTriple === 0 && nSingle === 0) v = 23.85; // shouldn't occur w/o H
        else if (h === 1 && nSingle === 2) v = in3Ring ? 21.94 : 12.03;
        else if (h === 1 && nDouble === 1) v = 23.85;
        else if (h === 2 && nSingle === 1) v = 26.02;
        else if (h === 3) v = 26.02;
        else v = 3.24;
      } else if (q === 1) {
        if (h === 0 && nSingle === 4) v = 0;
        else if (h === 0 && nSingle === 2 && nDouble === 1) v = 3.01;
        else if (h === 0 && nSingle === 1 && nTriple === 1) v = 4.36;
        else if (h === 1 && nSingle === 3) v = 4.44;
        else if (h === 1 && nDouble === 1) v = 13.97;
        else if (h === 2 && nSingle === 2) v = 16.61;
        else if (h === 2 && nDouble === 1) v = 25.59;
        else if (h === 3) v = 27.64;
        else v = 0;
      } else if (q === -1) v = 25.59;
    } else {
      if (isArom) v = 13.14;
      else if (q === 0) {
        if (h === 0 && nSingle === 2) v = in3Ring ? 12.53 : 9.23;
        else if (h === 0 && nDouble === 1) v = 17.07;
        else if (h >= 1) v = 20.23;
      } else if (q === -1) v = 23.06;
      else if (q === 1) v = h ? 0 : 0;
    }
    total += v;
  }
  return total;
}

/** Everything the properties panel needs for a set of atoms (default: whole molecule). */
export function computeProperties(mol, atomIds = null) {
  if (mol.isEmpty()) return null;
  const ids = atomIds ?? [...mol.atoms.keys()];
  const ringInfo = analyzeRings(mol);
  const counts = elementCounts(mol, ids);
  const charge = totalCharge(mol, ids);
  const mw = molecularWeight(counts);
  const hbd = hBondDonors(mol, ids);
  const hba = hBondAcceptors(mol, ids);
  const idSet = new Set(ids);
  const rings = ringInfo.rings.filter((r) => idSet.has(r[0]));
  const aromaticRings = ringInfo.aromatic.aromaticRings.filter((r) => idSet.has(r[0]));
  const rot = rotatableBonds(mol, ringInfo, ids);
  const psa = tpsa(mol, ringInfo, ids);
  const heavyAtoms = ids.filter((id) => mol.getAtom(id).el !== 'H').length;
  const lipinskiViolations = [mw > 500, hbd > 5, hba > 10].filter(Boolean).length;
  const exact = exactMass(counts, charge);
  return {
    counts, charge,
    formula: formulaString(counts, charge),
    formulaHTML: formulaHTML(counts, charge),
    molecularWeight: mw,
    exactMass: exact,
    mz: charge ? exact / Math.abs(charge) : exact,
    elementalAnalysis: elementalAnalysis(counts),
    isotopePattern: isotopePattern(counts, charge),
    degreeOfUnsaturation: degreeOfUnsaturation(counts),
    heavyAtoms,
    rings: rings.length,
    aromaticRings: aromaticRings.length,
    hBondDonors: hbd,
    hBondAcceptors: hba,
    rotatableBonds: rot,
    tpsa: psa,
    lipinskiViolations,
    fragments: mol.fragments().filter((f) => idSet.has(f[0])).length,
  };
}
