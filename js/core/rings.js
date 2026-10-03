// Ring perception (SSSR) and aromaticity detection.

import { valenceElectrons } from './elements.js';

/**
 * Smallest Set of Smallest Rings.
 * Returns an array of rings; each ring is an array of atom ids in cyclic order.
 */
export function findSSSR(mol) {
  const frags = mol.fragments();
  const nRings = mol.bondCount - mol.atomCount + frags.length;
  if (nRings <= 0) return [];

  // Candidate cycles: for every bond, the shortest path between its endpoints
  // that does not use the bond itself. Also add candidates from every pair of
  // shortest paths through each atom for better coverage of bridged systems.
  const candidates = new Map();
  const addCycle = (cycle) => {
    if (!cycle || cycle.length < 3) return;
    const key = [...cycle].sort((a, b) => a - b).join(',');
    if (!candidates.has(key)) candidates.set(key, cycle);
  };
  for (const bond of mol.bonds.values()) {
    const path = shortestPath(mol, bond.a1, bond.a2, bond.id);
    if (path) addCycle(path);
  }
  for (const atomId of mol.atoms.keys()) {
    if (mol.degree(atomId) < 3) continue;
    for (const c of cyclesThroughAtom(mol, atomId)) addCycle(c);
  }

  const sorted = [...candidates.values()].sort((a, b) => a.length - b.length);
  const bondIndex = new Map();
  let i = 0;
  for (const b of mol.bonds.values()) bondIndex.set(b.id, i++);

  // Gaussian elimination over GF(2) using bigint bitsets of bond indices.
  const basis = []; // {vec, pivot}
  const rings = [];
  for (const cycle of sorted) {
    let vec = cycleVector(mol, cycle, bondIndex);
    for (const { vec: bv, pivot } of basis) {
      if ((vec >> pivot) & 1n) vec ^= bv;
    }
    if (vec === 0n) continue;
    const pivot = highestBit(vec);
    basis.push({ vec, pivot });
    basis.sort((a, b) => (b.pivot > a.pivot ? 1 : -1));
    rings.push(cycle);
    if (rings.length >= nRings) break;
  }
  return rings;
}

function highestBit(v) {
  let p = 0n;
  while (v > 1n) { v >>= 1n; p++; }
  return p;
}

function cycleVector(mol, cycle, bondIndex) {
  let v = 0n;
  for (let k = 0; k < cycle.length; k++) {
    const b = mol.bondBetween(cycle[k], cycle[(k + 1) % cycle.length]);
    if (b) v |= 1n << BigInt(bondIndex.get(b.id));
  }
  return v;
}

function shortestPath(mol, from, to, excludeBondId) {
  const prev = new Map([[from, null]]);
  const queue = [from];
  while (queue.length) {
    const cur = queue.shift();
    if (cur === to) break;
    for (const b of mol.bondsOf(cur)) {
      if (b.id === excludeBondId) continue;
      const nb = b.a1 === cur ? b.a2 : b.a1;
      if (!prev.has(nb)) { prev.set(nb, cur); queue.push(nb); }
    }
  }
  if (!prev.has(to)) return null;
  const path = [];
  for (let c = to; c !== null; c = prev.get(c)) path.push(c);
  return path.reverse();
}

function cyclesThroughAtom(mol, root) {
  // BFS tree from root; every non-tree edge (u,v) where the paths root->u and
  // root->v only share root gives a cycle.
  const prev = new Map([[root, null]]);
  const depth = new Map([[root, 0]]);
  const queue = [root];
  const out = [];
  const pathTo = (n) => { const p = []; for (let c = n; c !== null; c = prev.get(c)) p.push(c); return p; };
  while (queue.length) {
    const cur = queue.shift();
    for (const nb of mol.neighbors(cur)) {
      if (!prev.has(nb)) {
        prev.set(nb, cur); depth.set(nb, depth.get(cur) + 1); queue.push(nb);
      } else if (prev.get(cur) !== nb && depth.get(nb) >= depth.get(cur)) {
        const pu = pathTo(cur), pv = pathTo(nb);
        const su = new Set(pu.slice(0, -1));
        if (pv.slice(0, -1).some((x) => su.has(x))) continue;
        out.push([...pu.reverse(), ...pv.slice(0, -1)]);
      }
    }
  }
  return out;
}

/** Map atomId -> array of ring indices; bondId -> array of ring indices. */
export function ringMembership(mol, rings) {
  const atomRings = new Map();
  const bondRings = new Map();
  rings.forEach((ring, idx) => {
    for (let k = 0; k < ring.length; k++) {
      const a = ring[k];
      if (!atomRings.has(a)) atomRings.set(a, []);
      atomRings.get(a).push(idx);
      const b = mol.bondBetween(a, ring[(k + 1) % ring.length]);
      if (b) {
        if (!bondRings.has(b.id)) bondRings.set(b.id, []);
        bondRings.get(b.id).push(idx);
      }
    }
  });
  return { atomRings, bondRings };
}

/**
 * Hückel aromaticity on SSSR rings (and fused pairs, e.g. azulene-like systems).
 * Returns { atoms: Set<id>, bonds: Set<id>, rings: number[][], aromaticRings: number[][] }.
 */
export function perceiveAromaticity(mol, rings = findSSSR(mol)) {
  const atoms = new Set();
  const bonds = new Set();
  const aromaticRings = [];

  const ringIsAromatic = (ring) => {
    if (ring.length < 3) return false;
    let pi = 0;
    for (let k = 0; k < ring.length; k++) {
      const c = piContribution(mol, ring[k], ring);
      if (c < 0) return false;
      pi += c;
    }
    return pi % 4 === 2;
  };

  const markRing = (ring) => {
    aromaticRings.push(ring);
    for (let k = 0; k < ring.length; k++) {
      atoms.add(ring[k]);
      const b = mol.bondBetween(ring[k], ring[(k + 1) % ring.length]);
      if (b) bonds.add(b.id);
    }
  };

  for (const ring of rings) if (ringIsAromatic(ring)) markRing(ring);

  // Fused bicyclic check (e.g. azulene: 5+7 rings each non-aromatic alone).
  for (let i = 0; i < rings.length; i++) {
    for (let j = i + 1; j < rings.length; j++) {
      const shared = rings[i].filter((a) => rings[j].includes(a));
      if (shared.length !== 2) continue;
      if (rings[i].every((a) => atoms.has(a)) && rings[j].every((a) => atoms.has(a))) continue;
      const envelope = fuseRings(rings[i], rings[j], shared);
      if (envelope && ringIsAromatic(envelope)) { markRing(rings[i]); markRing(rings[j]); }
    }
  }
  return { atoms, bonds, rings, aromaticRings };
}

function fuseRings(r1, r2, shared) {
  // Build the perimeter cycle of two rings sharing exactly one bond.
  const [s1, s2] = shared;
  const rot = (r, start) => { const i = r.indexOf(start); return [...r.slice(i), ...r.slice(0, i)]; };
  let a = rot(r1, s1);
  if (a[1] === s2) a = [a[0], ...a.slice(1).reverse()]; // ensure s2 is last
  if (a[a.length - 1] !== s2) return null;
  let b = rot(r2, s2);
  if (b[1] === s1) b = [b[0], ...b.slice(1).reverse()];
  if (b[b.length - 1] !== s1) return null;
  return [...a, ...b.slice(1, -1)];
}

/**
 * π-electron contribution of an atom to a ring, or -1 if the atom cannot be
 * part of an aromatic system (sp3 with no lone pair).
 */
function piContribution(mol, atomId, ring) {
  const atom = mol.getAtom(atomId);
  const bonds = mol.bondsOf(atomId);
  const ringSet = new Set(ring);
  let ringDouble = false;
  let exoDouble = null;
  for (const b of bonds) {
    if (b.order === 3) return -1;
    if (b.order === 2) {
      const other = b.a1 === atomId ? b.a2 : b.a1;
      if (ringSet.has(other)) ringDouble = true;
      else exoDouble = mol.getAtom(other);
    }
  }
  // Double bond to a neighbour belonging to another (fused) ring that is in
  // this ring set? Treat as ring double if both atoms are in ring.
  if (ringDouble) return 1;
  if (exoDouble) {
    // Exocyclic C=O, C=N, C=S etc. pull the electrons out of the ring.
    if (atom.el === 'C' && ['O', 'N', 'S'].includes(exoDouble.el)) return 0;
    if (atom.el === 'C') return 1; // exocyclic C=C (fulvene-like) - treat as 1
    return -1;
  }
  const ve = valenceElectrons(atom.el) - atom.charge;
  const totalBonds = mol.bondOrderSum(atomId) + mol.implicitH(atomId);
  const lonePairElectrons = ve - totalBonds - (atom.radical || 0);
  if (atom.el === 'C') {
    if (atom.charge === 1) return 0; // carbocation: empty p orbital
    if (atom.charge === -1) return 2; // carbanion
    if (atom.radical === 1) return 1;
    return -1; // sp3 carbon
  }
  if (atom.el === 'B') return atom.charge === -1 ? -1 : 0;
  if (lonePairElectrons >= 2) return 2; // pyrrole N, furan O, thiophene S
  if (lonePairElectrons === 1) return 1;
  return -1;
}

/** Convenience: SSSR + aromaticity + membership in one object. */
export function analyzeRings(mol) {
  const rings = findSSSR(mol);
  const membership = ringMembership(mol, rings);
  const aromatic = perceiveAromaticity(mol, rings);
  return { rings, ...membership, aromatic };
}
