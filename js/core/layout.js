// 2D coordinate generation (structure diagram layout).
//
// Strategy:
//  1. Perceive rings (SSSR) and group them into ring systems.
//  2. Lay out each ring system rigidly in local coordinates: the first ring
//     as a regular polygon, further rings fused / spiro / bridged onto it by
//     fitting the unplaced atoms of each ring on a circular arc.
//  3. Grow the molecule breadth-first from the largest ring system (or from
//     the end of the longest chain): chains zigzag at 120°, substituents are
//     spread in the largest free angular gap, triple bonds and cumulenes are
//     linear, ring systems are attached rigidly pointing away from the parent.
//  4. Relieve overlaps by reflecting / rotating substituents around acyclic
//     bonds, pick a pleasant orientation, and arrange fragments left-to-right.
//
// Coordinates use bond length 1.0 and a y-down axis (see molecule.js).

import { findSSSR } from './rings.js';

const TAU = Math.PI * 2;
const DEG = Math.PI / 180;
const FRAGMENT_GAP = 1.5;

/**
 * Assign 2D coordinates (bond length 1.0, y down) to every atom of `mol` in
 * place. Fragments are placed side by side and the drawing is centred on the
 * origin. Returns the molecule for chaining.
 * @param {import('./molecule.js').Molecule} mol
 */
export function layoutMolecule(mol) {
  if (mol.atomCount === 0) return mol;
  const rings = findSSSR(mol);
  const frags = mol.fragments();
  let x = 0;
  for (const frag of frags) {
    layoutAtoms(mol, frag, rings);
    const bb = mol.bbox(frag);
    mol.translate(x - bb.minX, -bb.cy, frag);
    x += bb.width + FRAGMENT_GAP;
  }
  const bb = mol.bbox();
  mol.translate(-bb.cx, -bb.cy);
  return mol;
}

/**
 * Lay out a single connected set of atoms in place, centred on the origin.
 * The atoms should form one or more complete fragments of `mol`.
 * @param {import('./molecule.js').Molecule} mol
 * @param {Iterable<number>} atomIds
 */
export function layoutFragment(mol, atomIds) {
  const ids = [...atomIds];
  if (!ids.length) return mol;
  const sub = mol.subMolecule(ids);
  layoutMolecule(sub);
  for (const id of ids) {
    const a = sub.getAtom(id);
    const t = mol.getAtom(id);
    t.x = a.x; t.y = a.y;
  }
  return mol;
}

// ------------------------------------------------------------------ helpers

const sub = (p, q) => ({ x: p.x - q.x, y: p.y - q.y });
const dist = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
const angleOf = (v) => Math.atan2(v.y, v.x);
const fromAngle = (a, len = 1) => ({ x: Math.cos(a) * len, y: Math.sin(a) * len });
const normAngle = (a) => ((a % TAU) + TAU) % TAU;

function rotateAround(p, c, ang) {
  const cos = Math.cos(ang), sin = Math.sin(ang);
  const dx = p.x - c.x, dy = p.y - c.y;
  return { x: c.x + dx * cos - dy * sin, y: c.y + dx * sin + dy * cos };
}

function reflectAcross(p, a, b) {
  const dx = b.x - a.x, dy = b.y - a.y;
  const len2 = dx * dx + dy * dy || 1;
  const t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
  const fx = a.x + t * dx, fy = a.y + t * dy;
  return { x: 2 * fx - p.x, y: 2 * fy - p.y };
}

/** Group rings sharing at least one atom into ring systems. */
function ringSystems(rings) {
  const parent = rings.map((_, i) => i);
  const find = (i) => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  const owner = new Map();
  rings.forEach((r, i) => {
    for (const a of r) {
      if (owner.has(a)) parent[find(i)] = find(owner.get(a));
      else owner.set(a, i);
    }
  });
  const groups = new Map();
  rings.forEach((r, i) => {
    const g = find(i);
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push(r);
  });
  return [...groups.values()].map((rs) => ({ rings: rs, atoms: new Set(rs.flat()) }));
}

// -------------------------------------------------------- ring system layout

function placeRingSystem(mol, rings) {
  const pos = new Map();
  const sets = rings.map((r) => new Set(r));
  const shared = (i, j) => rings[j].reduce((c, a) => c + (sets[i].has(a) ? 1 : 0), 0);
  // First ring: most fused neighbours, then size closest to 6.
  let first = 0, bestScore = -Infinity;
  rings.forEach((r, i) => {
    let fused = 0;
    for (let j = 0; j < rings.length; j++) if (j !== i && shared(i, j) >= 2) fused++;
    const score = fused * 10 - Math.abs(r.length - 6) - (r.length > 8 ? 5 : 0);
    if (score > bestScore) { bestScore = score; first = i; }
  });
  placeFirstPolygon(rings[first], pos);
  const done = new Set([first]);
  while (done.size < rings.length) {
    let best = -1, bestKey = -Infinity;
    rings.forEach((r, i) => {
      if (done.has(i)) return;
      const cnt = r.reduce((c, a) => c + (pos.has(a) ? 1 : 0), 0);
      // Prefer rings with exactly two placed atoms (simple fusion), then more.
      const key = (cnt === 2 ? 1000 : cnt * 10) - Math.abs(r.length - 6) * 0.1 - (cnt === 0 ? 1e6 : 0);
      if (key > bestKey) { bestKey = key; best = i; }
    });
    placeRing(mol, rings[best], pos);
    done.add(best);
  }
  if (needsRefinement(mol, pos)) stressRefine(mol, pos);
  return pos;
}

/** True when a ring system layout has distorted bonds or clashing atoms (cages, bridges). */
function needsRefinement(mol, pos) {
  const ids = [...pos.keys()];
  for (const id of ids) {
    for (const nb of mol.neighbors(id)) {
      if (pos.has(nb) && Math.abs(dist(pos.get(id), pos.get(nb)) - 1) > 0.05) return true;
    }
  }
  for (let i = 0; i < ids.length; i++) {
    for (let j = i + 1; j < ids.length; j++) {
      if (dist(pos.get(ids[i]), pos.get(ids[j])) < 0.75 && !mol.bondBetween(ids[i], ids[j])) return true;
    }
  }
  return false;
}

/** Stress majorization on topological distances (used for bridged / cage systems). */
function stressRefine(mol, pos) {
  const ids = [...pos.keys()];
  const n = ids.length;
  const index = new Map(ids.map((id, k) => [id, k]));
  const D = ids.map((id) => {
    const row = new Array(n).fill(Infinity);
    row[index.get(id)] = 0;
    const q = [id];
    while (q.length) {
      const c = q.shift();
      for (const nb of mol.neighbors(c)) {
        const k = index.get(nb);
        if (k === undefined || row[k] !== Infinity) continue;
        row[k] = row[index.get(c)] + 1;
        q.push(nb);
      }
    }
    return row.map((k) => (k === 0 ? 0 : k === 1 ? 1 : 1 + (k - 1) * 0.75));
  });
  const X = ids.map((id, k) => {
    const p = pos.get(id);
    // Deterministic jitter separates coincident atoms.
    return { x: p.x + Math.cos(k * 2.4) * 0.05, y: p.y + Math.sin(k * 2.4) * 0.05 };
  });
  for (let it = 0; it < 300; it++) {
    for (let i = 0; i < n; i++) {
      let sx = 0, sy = 0, sw = 0;
      for (let j = 0; j < n; j++) {
        if (i === j || !Number.isFinite(D[i][j])) continue;
        const t = D[i][j];
        const w = 1 / (t * t);
        const dx = X[i].x - X[j].x, dy = X[i].y - X[j].y;
        const d = Math.hypot(dx, dy) || 1e-6;
        sx += w * (X[j].x + (t * dx) / d);
        sy += w * (X[j].y + (t * dy) / d);
        sw += w;
      }
      if (sw) { X[i].x = sx / sw; X[i].y = sy / sw; }
    }
  }
  // Normalise to mean bond length 1.
  let sum = 0, cnt = 0;
  for (let i = 0; i < n; i++) {
    for (const nb of mol.neighbors(ids[i])) {
      const j = index.get(nb);
      if (j !== undefined && j > i) { sum += Math.hypot(X[i].x - X[j].x, X[i].y - X[j].y); cnt++; }
    }
  }
  const f = cnt ? cnt / sum : 1;
  ids.forEach((id, k) => pos.set(id, { x: X[k].x * f, y: X[k].y * f }));
}

function placeFirstPolygon(ring, pos) {
  const n = ring.length;
  const R = 1 / (2 * Math.sin(Math.PI / n));
  const start = n % 4 === 0 ? -Math.PI / 2 + Math.PI / n : -Math.PI / 2;
  ring.forEach((a, k) => pos.set(a, fromAngle(start + (k * TAU) / n, R)));
}

function placeRing(mol, ring, pos) {
  const n = ring.length;
  const placed = ring.map((a) => pos.has(a));
  const cnt = placed.filter(Boolean).length;
  if (cnt === n) return;
  if (cnt === 0) {
    // Disconnected from what is placed (should not happen): put it aside.
    const bb = boundsOf(pos);
    const R = 1 / (2 * Math.sin(Math.PI / n));
    ring.forEach((a, k) => {
      const p = fromAngle(-Math.PI / 2 + (k * TAU) / n, R);
      pos.set(a, { x: p.x + bb.maxX + R + 1, y: p.y });
    });
    return;
  }
  if (cnt === 1) {
    placeSpiro(mol, ring, pos);
    return;
  }
  const runs = [];
  for (let k = 0; k < n; k++) {
    if (placed[k] && !placed[(k + 1) % n]) {
      const run = [];
      let j = (k + 1) % n;
      while (!placed[j]) { run.push(ring[j]); j = (j + 1) % n; }
      runs.push({ p: ring[k], q: ring[j], run });
    }
  }
  for (const { p, q, run } of runs) placeArc(pos, p, q, run);
}

function placeSpiro(mol, ring, pos) {
  const n = ring.length;
  const si = ring.findIndex((a) => pos.has(a));
  const s = ring[si];
  const sp = pos.get(s);
  let dx = 0, dy = 0;
  for (const nb of mol.neighbors(s)) {
    const q = pos.get(nb);
    if (!q) continue;
    const d = dist(q, sp) || 1;
    dx -= (q.x - sp.x) / d; dy -= (q.y - sp.y) / d;
  }
  let len = Math.hypot(dx, dy);
  if (len < 1e-6) {
    const c = centroid([...pos.values()]);
    dx = sp.x - c.x; dy = sp.y - c.y; len = Math.hypot(dx, dy);
    if (len < 1e-6) { dx = 1; dy = 0; len = 1; }
  }
  dx /= len; dy /= len;
  const R = 1 / (2 * Math.sin(Math.PI / n));
  const c = { x: sp.x + dx * R, y: sp.y + dy * R };
  const a0 = angleOf(sub(sp, c));
  for (let k = 1; k < n; k++) {
    const a = ring[(si + k) % n];
    const p = fromAngle(a0 + (k * TAU) / n, R);
    pos.set(a, { x: c.x + p.x, y: c.y + p.y });
  }
}

/**
 * Place `run` (atoms between placed atoms p and q, in order) on a circular arc
 * with unit chords, bulging towards the less crowded side of the p-q line.
 */
function placeArc(pos, pId, qId, run) {
  const P = pos.get(pId), Q = pos.get(qId);
  const m = run.length + 1;
  const d = dist(P, Q);
  const M = { x: (P.x + Q.x) / 2, y: (P.y + Q.y) / 2 };
  const ux = d > 1e-9 ? (Q.x - P.x) / d : 1, uy = d > 1e-9 ? (Q.y - P.y) / d : 0;
  let nx = -uy, ny = ux;
  // Crowding on each side of the chord.
  let s1 = 0, s2 = 0, cx = 0, cy = 0, cn = 0;
  const R0 = Math.max(2.5, d * 1.2);
  for (const [id, a] of pos) {
    if (id === pId || id === qId) continue;
    cx += a.x; cy += a.y; cn++;
    const rx = a.x - M.x, ry = a.y - M.y;
    const dd = Math.hypot(rx, ry);
    if (dd > R0) continue;
    const side = rx * nx + ry * ny;
    const w = 1 / (dd + 0.3);
    if (side > 1e-3) s1 += w; else if (side < -1e-3) s2 += w;
  }
  let dir = 1;
  if (Math.abs(s1 - s2) > 1e-6) dir = s1 < s2 ? 1 : -1;
  else if (cn) dir = ((cx / cn - M.x) * nx + (cy / cn - M.y) * ny) > 0 ? -1 : 1;
  nx *= dir; ny *= dir;

  if (m <= d + 1e-6) {
    // Cannot span with unit chords: spread evenly on a slightly bowed line.
    run.forEach((a, i) => {
      const t = (i + 1) / m;
      const bow = Math.sin(Math.PI * t) * 0.3;
      pos.set(a, { x: P.x + (Q.x - P.x) * t + nx * bow, y: P.y + (Q.y - P.y) * t + ny * bow });
    });
    return;
  }
  const chordAng = (r) => 2 * Math.asin(Math.min(1, 1 / (2 * r)));
  const pqAng = (r) => 2 * Math.asin(Math.min(1, d / (2 * r)));
  const rMin = Math.max(d / 2, 0.5);
  const big = m * chordAng(rMin) + pqAng(rMin) >= TAU - 1e-9;
  const g = big ? (r) => m * chordAng(r) + pqAng(r) - TAU : (r) => m * chordAng(r) - pqAng(r);
  // g(rMin) >= 0 for big (decreasing), < 0 for small (increasing).
  let lo = rMin, hi = rMin;
  let guard = 0;
  while (guard++ < 200) {
    hi *= 2;
    if (big ? g(hi) < 0 : g(hi) > 0) break;
  }
  for (let it = 0; it < 80; it++) {
    const mid = (lo + hi) / 2;
    const v = g(mid);
    if (big ? v > 0 : v < 0) lo = mid; else hi = mid;
  }
  const r = (lo + hi) / 2;
  const h = Math.sqrt(Math.max(0, r * r - (d * d) / 4));
  const C = big ? { x: M.x + nx * h, y: M.y + ny * h } : { x: M.x - nx * h, y: M.y - ny * h };
  const th = chordAng(r);
  const aP = angleOf(sub(P, C));
  let bestSign = 1, bestErr = Infinity;
  for (const sign of [1, -1]) {
    const end = { x: C.x + r * Math.cos(aP + sign * th * m), y: C.y + r * Math.sin(aP + sign * th * m) };
    const mid = { x: C.x + r * Math.cos(aP + sign * th * m / 2), y: C.y + r * Math.sin(aP + sign * th * m / 2) };
    // Must end at Q and pass through the bulge side.
    const err = dist(end, Q) + (((mid.x - M.x) * nx + (mid.y - M.y) * ny) < 0 ? 10 : 0);
    if (err < bestErr) { bestErr = err; bestSign = sign; }
  }
  run.forEach((a, i) => {
    const ang = aP + bestSign * th * (i + 1);
    pos.set(a, { x: C.x + r * Math.cos(ang), y: C.y + r * Math.sin(ang) });
  });
}

function centroid(pts) {
  let x = 0, y = 0;
  for (const p of pts) { x += p.x; y += p.y; }
  return pts.length ? { x: x / pts.length, y: y / pts.length } : { x: 0, y: 0 };
}

function boundsOf(pos) {
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const p of pos.values()) {
    minX = Math.min(minX, p.x); maxX = Math.max(maxX, p.x);
    minY = Math.min(minY, p.y); maxY = Math.max(maxY, p.y);
  }
  return { minX, maxX, minY, maxY };
}

// ------------------------------------------------------------ fragment layout

function layoutAtoms(mol, ids, allRings) {
  const idSet = new Set(ids);
  if (ids.length === 1) {
    const a = mol.getAtom(ids[0]);
    a.x = 0; a.y = 0;
    return;
  }
  const rings = allRings.filter((r) => idSet.has(r[0]));
  const systems = ringSystems(rings);
  const sysOf = new Map();
  systems.forEach((s, i) => { for (const a of s.atoms) sysOf.set(a, i); });
  const ringBonds = new Set();
  const ringAngles = new Set(); // "n1|a|n2" with n1 < n2: angle inside a ring
  for (const r of rings) {
    for (let k = 0; k < r.length; k++) {
      const a = r[k], b = r[(k + 1) % r.length], c = r[(k + r.length - 1) % r.length];
      const bond = mol.bondBetween(a, b);
      if (bond) ringBonds.add(bond.id);
      ringAngles.add(angleKey(c, a, b));
    }
  }
  const local = systems.map((s) => placeRingSystem(mol, s.rings));
  const pos = new Map();
  const queue = [];
  const ctx = { mol, pos, ringAngles, idSet };

  if (systems.length) {
    let si = 0;
    systems.forEach((s, i) => { if (s.atoms.size > systems[si].atoms.size) si = i; });
    for (const [id, p] of local[si]) { pos.set(id, { ...p }); queue.push(id); }
  } else {
    const start = farthestAtom(mol, farthestAtom(mol, ids[0]));
    pos.set(start, { x: 0, y: 0 });
    const nbs = mol.neighbors(start);
    if (nbs.length === 1) {
      pos.set(nbs[0], fromAngle(-30 * DEG));
      queue.push(nbs[0]);
    } else {
      queue.push(start);
    }
  }

  const parent = new Map(); // chain atom -> atom it was placed from
  const sizeCache = new Map();
  const branchSize = (from, to) => {
    const key = `${from}>${to}`;
    if (!sizeCache.has(key)) {
      const seen = new Set([from, to]);
      const st = [to];
      while (st.length) {
        const c = st.pop();
        for (const nb of mol.neighbors(c)) if (!seen.has(nb)) { seen.add(nb); st.push(nb); }
      }
      sizeCache.set(key, seen.size - 1);
    }
    return sizeCache.get(key);
  };

  while (queue.length) {
    const A = queue.shift();
    const nbs = mol.neighbors(A);
    const U = nbs.filter((x) => !pos.has(x));
    if (!U.length) continue;
    const placedN = nbs.filter((x) => pos.has(x));
    const a = pos.get(A);
    let targets; // array of positions, one per slot
    let order = U;

    if (placedN.length === 1 && U.length === 1) {
      const P = placedN[0];
      const p = pos.get(P);
      const dv = sub(a, p);
      const dl = Math.hypot(dv.x, dv.y) || 1;
      const din = angleOf(dv);
      if (isLinear(mol, A)) {
        targets = [{ x: a.x + dv.x / dl, y: a.y + dv.y / dl }];
      } else {
        const c1 = { x: a.x + Math.cos(din + 60 * DEG), y: a.y + Math.sin(din + 60 * DEG) };
        const c2 = { x: a.x + Math.cos(din - 60 * DEG), y: a.y + Math.sin(din - 60 * DEG) };
        // Zigzag reference: the atom P was grown from, else P's other neighbours.
        const refs = parent.has(P) && parent.get(P) !== A ? [pos.get(parent.get(P))]
          : mol.neighbors(P).filter((x) => x !== A && pos.has(x)).map((x) => pos.get(x));
        let pref = c1, other = c2;
        if (refs.length) {
          const rc = centroid(refs);
          const sideRef = cross(dv, sub(rc, p));
          const side1 = cross(dv, sub(c1, p));
          if (Math.abs(sideRef) > 1e-6 && Math.sign(side1) === Math.sign(sideRef)) { pref = c2; other = c1; }
        }
        if (crowd(pos, pref, A) > crowd(pos, other, A) + 0.3) [pref, other] = [other, pref];
        targets = [pref];
      }
    } else {
      const k = U.length;
      const slots = [];
      if (placedN.length === 0) {
        for (let i = 0; i < k; i++) slots.push(-30 * DEG + (i * TAU) / k);
      } else if (placedN.length === 1) {
        const back = angleOf(sub(pos.get(placedN[0]), a));
        for (let i = 0; i < k; i++) slots.push(back + (TAU * (i + 1)) / (k + 1));
      } else {
        const { start, size } = freeGap(ctx, A, placedN);
        for (let i = 0; i < k; i++) slots.push(start + (size * (i + 1)) / (k + 1));
      }
      const slotPts = slots.map((ang) => ({ x: a.x + Math.cos(ang), y: a.y + Math.sin(ang) }));
      // Least crowded slot gets the biggest branch.
      const slotOrder = slotPts.map((p, i) => ({ i, c: crowd(pos, p, A) })).sort((x, y) => x.c - y.c || x.i - y.i);
      order = [...U].sort((x, y) => branchSize(A, y) - branchSize(A, x) || x - y);
      targets = slotOrder.map((s) => slotPts[s.i]);
    }

    order.forEach((u, i) => {
      const t = targets[i];
      if (sysOf.has(u)) {
        const si = sysOf.get(u);
        attachSystem(ctx, local[si], u, t, A);
        for (const id of local[si].keys()) queue.push(id);
      } else {
        pos.set(u, t);
        parent.set(u, A);
        queue.push(u);
      }
    });
  }

  // Any atom not reached (should not happen) - park at origin.
  for (const id of ids) if (!pos.has(id)) pos.set(id, { x: 0, y: 0 });

  relieveOverlaps(mol, ids, pos, ringBonds);
  orient(pos);
  for (const [id, p] of pos) {
    const at = mol.getAtom(id);
    at.x = p.x; at.y = p.y;
  }
}

function angleKey(n1, a, n2) {
  return n1 < n2 ? `${n1}|${a}|${n2}` : `${n2}|${a}|${n1}`;
}

const cross = (u, v) => u.x * v.y - u.y * v.x;

/** Linear geometry: triple bond, or two cumulated double bonds, with degree 2. */
function isLinear(mol, id) {
  const bonds = mol.bondsOf(id);
  if (bonds.length !== 2) return false;
  let doubles = 0;
  for (const b of bonds) {
    if (b.order >= 3) return true;
    if (b.order === 2) doubles++;
  }
  return doubles === 2;
}

/** Crowding score of a candidate point against already placed atoms. */
function crowd(pos, p, exclude) {
  let s = 0;
  for (const [id, q] of pos) {
    if (id === exclude) continue;
    const d = dist(p, q);
    if (d < 2.5) s += 1 / (d * d + 0.05);
  }
  return s;
}

/**
 * Largest free angular gap around atom `A` given placed neighbours, preferring
 * gaps that are not ring interiors. Returns { start, size } in radians.
 */
function freeGap(ctx, A, placedN) {
  const { pos, ringAngles } = ctx;
  const a = pos.get(A);
  const items = placedN.map((id) => ({ id, ang: normAngle(angleOf(sub(pos.get(id), a))) }))
    .sort((x, y) => x.ang - y.ang);
  let best = null;
  for (let i = 0; i < items.length; i++) {
    const cur = items[i], next = items[(i + 1) % items.length];
    let size = next.ang - cur.ang;
    if (i === items.length - 1) size += TAU;
    if (items.length === 1) size = TAU;
    const inRing = items.length > 1 && ringAngles.has(angleKey(cur.id, A, next.id));
    const score = size - (inRing ? 10 : 0);
    if (!best || score > best.score + 1e-9) best = { start: cur.ang, size, score };
  }
  return best;
}

/** Rigidly place a ring system so that atom u sits at `target` facing atom A. */
function attachSystem(ctx, loc, u, target, A) {
  const { mol, pos, ringAngles } = ctx;
  const lu = loc.get(u);
  const ringNbs = mol.neighbors(u).filter((x) => loc.has(x));
  const exo = mol.neighbors(u).filter((x) => !loc.has(x)).length;
  // Free gap in local coordinates.
  const tmp = { pos: loc, ringAngles };
  const { start, size } = freeGap(tmp, u, ringNbs);
  const slot = start + size / (exo + 1);
  const want = angleOf(sub(pos.get(A), target));
  const rot = want - slot;
  const axisB = { x: lu.x + Math.cos(slot), y: lu.y + Math.sin(slot) };
  let best = null;
  for (const mirror of [false, true]) {
    const out = new Map();
    let pen = 0;
    for (const [id, p0] of loc) {
      let p = mirror ? reflectAcross(p0, lu, axisB) : p0;
      p = rotateAround(p, lu, rot);
      const q = { x: p.x - lu.x + target.x, y: p.y - lu.y + target.y };
      out.set(id, q);
      for (const [pid, pp] of pos) {
        if (pid === A && id === u) continue;
        const d = dist(q, pp);
        if (d < 1.6) pen += (1.6 - d) ** 2;
      }
    }
    if (!best || pen < best.pen - 1e-9) best = { pen, out };
  }
  for (const [id, q] of best.out) pos.set(id, q);
}

function farthestAtom(mol, start) {
  const seen = new Map([[start, 0]]);
  const q = [start];
  let last = start;
  while (q.length) {
    const c = q.shift();
    last = c;
    for (const nb of mol.neighbors(c)) if (!seen.has(nb)) { seen.set(nb, seen.get(c) + 1); q.push(nb); }
  }
  return last;
}

// ---------------------------------------------------------- overlap relief

function relieveOverlaps(mol, ids, pos, ringBonds) {
  if (ids.length < 4) return;
  const LIMIT = 1.0;
  const idSet = new Set(ids);
  const nbrs = new Map(ids.map((id) => [id, new Set(mol.neighbors(id))]));
  const bonds = mol.bondList().filter((b) => idSet.has(b.a1) && !ringBonds.has(b.id));
  const pairPen = (i, pi, j, pj) => {
    const dx = pi.x - pj.x, dy = pi.y - pj.y;
    if (dx > LIMIT || dx < -LIMIT || dy > LIMIT || dy < -LIMIT) return 0;
    const d = Math.sqrt(dx * dx + dy * dy);
    if (d >= LIMIT || nbrs.get(i).has(j)) return 0;
    return (LIMIT - d) ** 2 + (d < 0.5 ? 1 : 0);
  };
  const clashPairs = () => {
    const out = [];
    for (let i = 0; i < ids.length; i++) {
      const pi = pos.get(ids[i]);
      for (let j = i + 1; j < ids.length; j++) {
        if (pairPen(ids[i], pi, ids[j], pos.get(ids[j])) > 0) out.push([ids[i], ids[j]]);
      }
    }
    return out;
  };
  const sideOf = (from, excludeTo) => {
    const seen = new Set([from]);
    const st = [from];
    while (st.length) {
      const c = st.pop();
      for (const nb of nbrs.get(c)) {
        if (c === from && nb === excludeTo) continue;
        if (!seen.has(nb)) { seen.add(nb); st.push(nb); }
      }
    }
    return seen;
  };
  const crossPen = (side, others, getPos) => {
    let s = 0;
    for (const i of side) {
      const pi = getPos(i);
      for (const j of others) s += pairPen(i, pi, j, pos.get(j));
    }
    return s;
  };
  for (let pass = 0; pass < 50; pass++) {
    const clashes = clashPairs();
    if (!clashes.length) return;
    let improved = false;
    for (const b of bonds) {
      let side = sideOf(b.a2, b.a1);
      let pivot = b.a1, head = b.a2;
      if (side.has(b.a1)) continue; // not a bridge after all
      if (side.size > ids.length / 2) { side = sideOf(b.a1, b.a2); pivot = b.a2; head = b.a1; }
      // Only bonds separating a clashing pair can help.
      if (!clashes.some(([x, y]) => side.has(x) !== side.has(y))) continue;
      const others = ids.filter((j) => !side.has(j));
      const e0 = crossPen(side, others, (i) => pos.get(i));
      if (e0 === 0) continue;
      const pp = pos.get(pivot), hp = pos.get(head);
      const candidates = [{ f: (p) => reflectAcross(p, pp, hp), cost: 0 }];
      for (const ang of [30, -30, 60, -60, 90, -90, 120, -120]) {
        candidates.push({ f: (p) => rotateAround(p, pp, ang * DEG), cost: Math.abs(ang) * 0.0005 });
        candidates.push({ f: (p) => rotateAround(reflectAcross(p, pp, hp), pp, ang * DEG), cost: Math.abs(ang) * 0.0005 + 0.0001 });
      }
      let best = null;
      for (const c of candidates) {
        const moved = new Map();
        for (const i of side) moved.set(i, c.f(pos.get(i)));
        const e = crossPen(side, others, (i) => moved.get(i)) + c.cost;
        if (e < e0 - 1e-6 && (!best || e < best.e)) best = { e, moved };
      }
      if (best) {
        for (const [i, p] of best.moved) pos.set(i, p);
        improved = true;
        break; // clash list is stale; recompute
      }
    }
    if (!improved) break;
  }
}

/** Rotate by a multiple of 30° when that makes the drawing clearly wider than tall. */
function orient(pos) {
  const pts = [...pos.values()];
  if (pts.length < 3) return;
  const c = centroid(pts);
  const score = (ang) => {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of pts) {
      const q = rotateAround(p, c, ang);
      minX = Math.min(minX, q.x); maxX = Math.max(maxX, q.x);
      minY = Math.min(minY, q.y); maxY = Math.max(maxY, q.y);
    }
    return (maxX - minX) - (maxY - minY);
  };
  let bestAng = 0, best = score(0);
  for (let k = 1; k < 6; k++) {
    const s = score(k * 30 * DEG);
    if (s > best + 0.75) { best = s; bestAng = k * 30 * DEG; }
  }
  if (bestAng) for (const [id, p] of pos) pos.set(id, rotateAround(p, c, bestAng));
}
