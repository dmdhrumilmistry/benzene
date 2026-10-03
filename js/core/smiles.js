// SMILES reader / writer.
//
// The parser supports the practical OpenSMILES subset: organic-subset and
// bracket atoms (isotope, chirality, H count, charge, atom class), all bond
// symbols, branches, ring closures (including %nn) and dot-separated
// fragments. Aromatic input is kekulized with a maximum matching.

import { Molecule, adjustValence, computeImplicitH } from './molecule.js';
import { isElement, VALENCES, ELEMENTS } from './elements.js';
import { perceiveAromaticity } from './rings.js';
import { layoutMolecule } from './layout.js';

const ORGANIC = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']);
const AROMATIC_ORGANIC = new Set(['b', 'c', 'n', 'o', 'p', 's']);
const AROMATIC_BRACKET = new Set(['b', 'c', 'n', 'o', 'p', 's', 'se', 'as', 'te', 'si', 'ge', 'sb', 'bi']);
const AROMATIC_WRITABLE = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'Se', 'As', 'Te']);
const SMILES_VALENCES = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
const BOND_CHARS = { '-': 1, '=': 2, '#': 3, $: 4, ':': 'ar', '/': 1, '\\': 1 };

// ------------------------------------------------------------------- parser

/**
 * Parse a SMILES string into a Molecule with 2D coordinates.
 * Anything after the first whitespace (e.g. a name) is ignored.
 * @param {string} str
 * @returns {Molecule}
 * @throws {Error} on invalid SMILES
 */
export function parseSmiles(str) {
  const mol = new Molecule();
  const s = String(str ?? '').trim().split(/\s+/)[0] || '';
  if (!s) return mol;

  const aromaticAtoms = new Set();
  const aromaticBonds = new Set();
  const ringOpen = new Map(); // ring number -> { atom, bond, pos }
  const stack = [];
  let prev = null;
  let pendingBond = null; // { order, ch, pos }
  let i = 0;

  const fail = (msg, at = i) => {
    throw new Error(`Invalid SMILES: ${msg} at position ${at + 1}\n  ${s}\n  ${' '.repeat(at)}^`);
  };

  const connect = (a, b, bondSpec, at) => {
    if (a === b) fail('atom bonded to itself', at);
    if (mol.bondBetween(a, b)) fail('duplicate bond between the same atoms', at);
    let order = 1;
    let arom = false;
    if (bondSpec) {
      if (bondSpec.order === 'ar') arom = true;
      else order = bondSpec.order;
    } else if (aromaticAtoms.has(a) && aromaticAtoms.has(b)) {
      arom = true;
    }
    const bond = mol.addBond(a, b, order);
    if (arom) aromaticBonds.add(bond.id);
    return bond;
  };

  const addAtom = (props, aromatic) => {
    const atom = mol.addAtom(props.atom);
    if (props.chirality) atom.chirality = props.chirality;
    if (aromatic) aromaticAtoms.add(atom.id);
    if (prev !== null) connect(prev, atom.id, pendingBond, pendingBond ? pendingBond.pos : i);
    else if (pendingBond) fail('bond without a preceding atom', pendingBond.pos);
    pendingBond = null;
    prev = atom.id;
  };

  while (i < s.length) {
    const ch = s[i];
    if (ch === '[') {
      const end = s.indexOf(']', i);
      if (end < 0) fail('unclosed bracket atom');
      const parsed = parseBracket(s.slice(i + 1, end), i, fail);
      addAtom(parsed, parsed.aromatic);
      i = end + 1;
    } else if (/[A-Za-z*]/.test(ch)) {
      const two = s.slice(i, i + 2);
      let sym;
      if (two === 'Cl' || two === 'Br') sym = two;
      else sym = ch;
      if (sym === '*') fail('wildcard atoms (*) are not supported');
      if (ORGANIC.has(sym)) {
        addAtom({ atom: { el: sym } }, false);
      } else if (AROMATIC_ORGANIC.has(sym)) {
        addAtom({ atom: { el: sym.toUpperCase() } }, true);
      } else {
        fail(`'${sym}' is not an organic-subset atom (use brackets, e.g. [${sym}])`);
      }
      i += sym.length;
    } else if (ch in BOND_CHARS) {
      if (pendingBond) fail('two consecutive bond symbols');
      if (prev === null) fail('bond without a preceding atom');
      pendingBond = { order: BOND_CHARS[ch], ch, pos: i };
      i++;
    } else if (ch === '(') {
      if (prev === null) fail('branch without a preceding atom');
      if (pendingBond) fail('bond symbol before branch');
      stack.push(prev);
      i++;
      if (s[i] === ')') fail('empty branch');
    } else if (ch === ')') {
      if (!stack.length) fail('unbalanced parenthesis');
      if (pendingBond) fail('bond symbol at end of branch');
      prev = stack.pop();
      i++;
    } else if (/[0-9%]/.test(ch)) {
      if (prev === null) fail('ring closure without a preceding atom');
      let num;
      const at = i;
      if (ch === '%') {
        if (s[i + 1] === '(') {
          const close = s.indexOf(')', i);
          if (close < 0) fail('unclosed %( ring number');
          num = Number(s.slice(i + 2, close));
          i = close + 1;
        } else {
          const digits = s.slice(i + 1, i + 3);
          if (!/^\d\d$/.test(digits)) fail("'%' must be followed by two digits");
          num = Number(digits);
          i += 3;
        }
      } else {
        num = Number(ch);
        i++;
      }
      if (ringOpen.has(num)) {
        const open = ringOpen.get(num);
        ringOpen.delete(num);
        let spec = pendingBond || open.bond;
        if (pendingBond && open.bond && pendingBond.order !== open.bond.order
          && !(pendingBond.order === 1 && open.bond.order === 1)) {
          fail(`conflicting bond orders for ring closure ${num}`, at);
        }
        connect(open.atom, prev, spec, at);
      } else {
        ringOpen.set(num, { atom: prev, bond: pendingBond, pos: at });
      }
      pendingBond = null;
    } else if (ch === '.') {
      if (pendingBond) fail('bond symbol before dot');
      if (stack.length) fail('dot inside a branch');
      prev = null;
      i++;
    } else {
      fail(`unexpected character '${ch}'`);
    }
  }
  if (pendingBond) fail('SMILES ends with a bond symbol', pendingBond.pos);
  if (stack.length) fail('unclosed branch', s.length - 1);
  if (ringOpen.size) {
    const [num, open] = ringOpen.entries().next().value;
    fail(`unclosed ring ${num}`, open.pos);
  }

  if (aromaticBonds.size) {
    const res = kekulize(mol, aromaticBonds, { aromaticAtoms });
    if (!res.success) {
      throw new Error(`Invalid SMILES: cannot kekulize aromatic system (check for missing [nH]) in ${s}`);
    }
  }
  layoutMolecule(mol);
  return mol;
}

function parseBracket(body, offset, fail) {
  const m = /^(\d+)?([A-Z][a-z]?|se|as|te|si|ge|sb|bi|[bcnops]|\*)(@@|@(?:TH[12]|AL[12]|SP[1-3]|TB\d{1,2}|OH\d{1,2})?)?(H\d*)?((?:[+-]\d+)|\+{1,3}|-{1,3})?(?::(\d+))?$/.exec(body);
  if (!m) fail(`malformed bracket atom [${body}]`, offset);
  const [, iso, sym, chir, hs, chg, cls] = m;
  if (sym === '*') fail('wildcard atoms (*) are not supported', offset);
  let el = sym;
  let aromatic = false;
  if (sym[0] === sym[0].toLowerCase()) {
    if (!AROMATIC_BRACKET.has(sym)) fail(`unknown aromatic atom '${sym}'`, offset);
    aromatic = true;
    el = sym[0].toUpperCase() + sym.slice(1);
  }
  if (!isElement(el)) fail(`unknown element '${el}'`, offset);
  let charge = 0;
  if (chg) {
    if (/^[+-]\d+$/.test(chg)) charge = Number(chg);
    else charge = (chg[0] === '+' ? 1 : -1) * chg.length;
  }
  const hCount = hs ? (hs.length > 1 ? Number(hs.slice(1)) : 1) : 0;
  const atom = { el, charge, hCount, isotope: iso ? Number(iso) : null, mapNo: cls ? Number(cls) : 0 };
  const chirality = chir === '@' || chir === '@@' ? chir : null;
  return { atom, aromatic, chirality };
}

// ------------------------------------------------------------- kekulization

/**
 * Convert aromatic bonds into alternating single/double bonds in place.
 * Bonds in `aromaticBondIds` are treated as order 1 ("aromatic"); atoms that
 * still have a free valence receive exactly one double bond via a maximum
 * matching (Edmonds' blossom). Aromatic bonds that are not in a ring become
 * single. If a perfect matching fails, uncharged 2-connected N/P atoms with
 * implicit H are tried as pyrrole-type [nH].
 * @param {Molecule} mol
 * @param {Iterable<number>} aromaticBondIds
 * @param {{aromaticAtoms?: Set<number>}} [opts]
 * @returns {{success: boolean, unmatched: number[]}}
 */
export function kekulize(mol, aromaticBondIds, opts = {}) {
  const arom = new Set();
  for (const id of aromaticBondIds) {
    const b = mol.getBond(id);
    if (b) { b.order = 1; arom.add(id); }
  }
  // Aromatic bonds that are bridges (not in any ring) are plain single bonds.
  for (const id of bridges(mol)) arom.delete(id);
  const atoms = new Set();
  for (const id of arom) { const b = mol.getBond(id); atoms.add(b.a1); atoms.add(b.a2); }
  if (!atoms.size) return { success: true, unmatched: [] };

  const fixedH = new Set(); // atoms we turned into [nH]
  const needs = (id) => {
    const atom = mol.getAtom(id);
    let used = 0;
    for (const b of mol.bondsOf(id)) used += b.order;
    if (atom.hCount !== null && atom.hCount !== undefined) used += atom.hCount;
    used += atom.radical || 0;
    const vals = VALENCES[atom.el];
    if (!vals) return false;
    for (const v of vals) {
      const av = adjustValence(atom.el, v, atom.charge);
      if (av >= used) return av - used >= 1;
    }
    return false;
  };

  const solve = () => {
    const nodes = [...atoms].filter((a) => needs(a));
    const index = new Map(nodes.map((a, k) => [a, k]));
    const adj = nodes.map(() => []);
    const edgeBond = new Map();
    for (const id of arom) {
      const b = mol.getBond(id);
      const u = index.get(b.a1), v = index.get(b.a2);
      if (u === undefined || v === undefined) continue;
      adj[u].push(v); adj[v].push(u);
      edgeBond.set(`${Math.min(u, v)},${Math.max(u, v)}`, b);
    }
    const match = maxMatching(nodes.length, adj);
    const unmatched = nodes.filter((_, k) => match[k] < 0);
    return { nodes, match, edgeBond, unmatched };
  };

  let res = solve();
  if (res.unmatched.length) {
    // Try pyrrole-type nitrogens/phosphorus one at a time.
    const candidates = [...atoms].filter((id) => {
      const a = mol.getAtom(id);
      return (a.el === 'N' || a.el === 'P') && a.charge === 0 && a.hCount === null
        && mol.degree(id) === 2 && (!opts.aromaticAtoms || opts.aromaticAtoms.has(id));
    });
    let progress = true;
    while (res.unmatched.length && progress) {
      progress = false;
      for (const c of candidates) {
        if (fixedH.has(c)) continue;
        mol.getAtom(c).hCount = 1;
        const r2 = solve();
        if (r2.unmatched.length < res.unmatched.length) {
          fixedH.add(c); res = r2; progress = true;
          break;
        }
        mol.getAtom(c).hCount = null;
      }
    }
  }
  const { match, edgeBond } = res;
  for (let u = 0; u < match.length; u++) {
    const v = match[u];
    if (v > u) edgeBond.get(`${u},${v}`).order = 2;
  }
  mol.invalidate();
  return { success: res.unmatched.length === 0, unmatched: res.unmatched };
}

/** Bond ids that are bridges (removal disconnects the graph). */
function bridges(mol) {
  const disc = new Map();
  const low = new Map();
  const out = [];
  let t = 0;
  for (const root of mol.atoms.keys()) {
    if (disc.has(root)) continue;
    // Iterative DFS: stack of [atom, parentBondId, bondIterIndex]
    const stack = [[root, -1, 0]];
    disc.set(root, t); low.set(root, t); t++;
    while (stack.length) {
      const top = stack[stack.length - 1];
      const [v, pb] = top;
      const bonds = mol.bondsOf(v);
      if (top[2] < bonds.length) {
        const b = bonds[top[2]++];
        if (b.id === pb) continue;
        const w = b.a1 === v ? b.a2 : b.a1;
        if (disc.has(w)) {
          low.set(v, Math.min(low.get(v), disc.get(w)));
        } else {
          disc.set(w, t); low.set(w, t); t++;
          stack.push([w, b.id, 0]);
        }
      } else {
        stack.pop();
        if (stack.length) {
          const u = stack[stack.length - 1][0];
          low.set(u, Math.min(low.get(u), low.get(v)));
          if (low.get(v) > disc.get(u)) out.push(pb);
        }
      }
    }
  }
  return out;
}

/** Maximum matching in a general graph (Edmonds' blossom algorithm). */
function maxMatching(n, adj) {
  const match = new Array(n).fill(-1);
  const p = new Array(n);
  const base = new Array(n);
  const used = new Array(n);
  const blossom = new Array(n);

  const lca = (a, b) => {
    const seen = new Array(n).fill(false);
    for (;;) {
      a = base[a]; seen[a] = true;
      if (match[a] === -1) break;
      a = p[match[a]];
    }
    for (;;) {
      b = base[b];
      if (seen[b]) return b;
      b = p[match[b]];
    }
  };
  const markPath = (v, b, child) => {
    while (base[v] !== b) {
      blossom[base[v]] = blossom[base[match[v]]] = true;
      p[v] = child;
      child = match[v];
      v = p[match[v]];
    }
  };
  const findPath = (root) => {
    used.fill(false); p.fill(-1);
    for (let i = 0; i < n; i++) base[i] = i;
    used[root] = true;
    const q = [root];
    let qh = 0;
    while (qh < q.length) {
      const v = q[qh++];
      for (const to of adj[v]) {
        if (base[v] === base[to] || match[v] === to) continue;
        if (to === root || (match[to] !== -1 && p[match[to]] !== -1)) {
          const cur = lca(v, to);
          blossom.fill(false);
          markPath(v, cur, to);
          markPath(to, cur, v);
          for (let i = 0; i < n; i++) {
            if (blossom[base[i]]) {
              base[i] = cur;
              if (!used[i]) { used[i] = true; q.push(i); }
            }
          }
        } else if (p[to] === -1) {
          p[to] = v;
          if (match[to] === -1) return to;
          used[match[to]] = true;
          q.push(match[to]);
        }
      }
    }
    return -1;
  };

  // Greedy start: lowest-degree vertices first.
  const order = [...Array(n).keys()].sort((a, b) => adj[a].length - adj[b].length);
  for (const v of order) {
    if (match[v] !== -1) continue;
    const cands = adj[v].filter((w) => match[w] === -1).sort((a, b) => adj[a].length - adj[b].length);
    if (cands.length) { match[v] = cands[0]; match[cands[0]] = v; }
  }
  for (let v = 0; v < n; v++) {
    if (match[v] !== -1) continue;
    let u = findPath(v);
    while (u !== -1) {
      const pv = p[u];
      const ppv = match[pv];
      match[u] = pv; match[pv] = u;
      u = ppv;
    }
  }
  return match;
}

// ------------------------------------------------------------------- writer

/**
 * Write a (canonical-ish, deterministic) SMILES string.
 * @param {Molecule} mol
 * @param {{aromatic?: boolean}} [opts] aromatic: write aromatic rings in lowercase
 * @returns {string}
 */
export function writeSmiles(mol, { aromatic = true } = {}) {
  if (!mol || mol.atomCount === 0) return '';
  const aroAtoms = new Set();
  const aroBonds = new Set();
  if (aromatic) {
    const perc = perceiveAromaticity(mol);
    for (const a of perc.atoms) if (AROMATIC_WRITABLE.has(mol.getAtom(a).el)) aroAtoms.add(a);
    for (const b of perc.bonds) aroBonds.add(b);
    pruneAromatic(mol, aroAtoms, aroBonds);
  }
  const ranks = rankAtoms(mol, aroAtoms);
  const visited = new Set();
  const parts = [];
  const frags = mol.fragments().map((f) => f.sort((a, b) => ranks.get(a) - ranks.get(b)));
  // Larger fragments first, then by best rank.
  frags.sort((a, b) => b.length - a.length || ranks.get(a[0]) - ranks.get(b[0]));
  for (const frag of frags) {
    // Prefer a terminal atom; frag is sorted by rank so the first match wins.
    const start = frag.find((id) => mol.degree(id) === 1) ?? frag[0];
    parts.push(writeFragment(mol, start, ranks, aroAtoms, aroBonds, visited));
  }
  return parts.join('.');
}

/** Drop aromatic atoms whose lowercase form would not re-kekulize to the same structure. */
function pruneAromatic(mol, aroAtoms, aroBonds) {
  const isAroBond = (b) => aroBonds.has(b.id) && aroAtoms.has(b.a1) && aroAtoms.has(b.a2);
  for (let iter = 0; iter < 50; iter++) {
    for (const id of [...aroBonds]) {
      const b = mol.getBond(id);
      if (!b || !aroAtoms.has(b.a1) || !aroAtoms.has(b.a2)) aroBonds.delete(id);
    }
    let changed = false;
    for (const id of [...aroAtoms]) {
      const atom = mol.getAtom(id);
      const bonds = mol.bondsOf(id);
      const aroB = bonds.filter(isAroBond);
      if (aroB.length < 2) { aroAtoms.delete(id); changed = true; continue; }
      const hasDouble = aroB.some((b) => b.order === 2);
      if (aroB.some((b) => b.order > 2)) { aroAtoms.delete(id); changed = true; continue; }
      let used = aroB.length;
      for (const b of bonds) if (!isAroBond(b)) used += b.order;
      // Would the organic or bracket form request a double bond?
      const organic = canWriteOrganic(mol, id, true, used, hasDouble);
      const h = mol.implicitH(id);
      const need = needsDouble(atom.el, atom.charge, used + h + (atom.radical || 0));
      if (!organic && need !== hasDouble) { aroAtoms.delete(id); changed = true; }
    }
    if (!changed) break;
  }
}

function needsDouble(el, charge, used) {
  const vals = VALENCES[el];
  if (!vals) return false;
  for (const v of vals) {
    const av = adjustValence(el, v, charge);
    if (av >= used) return av - used >= 1;
  }
  return false;
}

/** Implicit H that the organic-subset form would imply. */
function smilesImplicitH(el, used) {
  const vals = SMILES_VALENCES[el];
  if (!vals) return 0;
  for (const v of vals) if (v >= used) return v - used;
  return 0;
}

/**
 * Can atom `id` be written without brackets? For aromatic atoms `used` is the
 * valence with aromatic bonds counted as 1 and `hasDouble` whether its Kekulé
 * form has a ring double bond.
 */
function canWriteOrganic(mol, id, aromatic, used, hasDouble) {
  const atom = mol.getAtom(id);
  if (!ORGANIC.has(atom.el) || atom.charge || atom.isotope || atom.radical || atom.mapNo) return false;
  const h = mol.implicitH(id);
  if (!aromatic) return smilesImplicitH(atom.el, mol.bondOrderSum(id)) === h;
  if (!AROMATIC_ORGANIC.has(atom.el.toLowerCase())) return false;
  // Parser: organic aromatic atom gets a double bond iff it has free valence,
  // then implicit H from the normal valence model.
  const need = needsDouble(atom.el, 0, used);
  if (need !== hasDouble) return false;
  return computeImplicitH(atom.el, 0, used + (hasDouble ? 1 : 0)) === h;
}

/** Morgan-like ranking: Map atomId -> rank (0 = first). */
function rankAtoms(mol, aroAtoms) {
  const ids = [...mol.atoms.keys()];
  const inv = new Map();
  for (const id of ids) {
    const a = mol.getAtom(id);
    const z = ELEMENTS[a.el]?.z ?? 0;
    inv.set(id, [mol.degree(id), z, mol.implicitH(id), a.charge + 8, a.isotope ?? 0, aroAtoms.has(id) ? 1 : 0]);
  }
  const toRanks = (keyOf) => {
    const sorted = [...ids].sort((x, y) => cmp(keyOf(x), keyOf(y)));
    const r = new Map();
    let rank = 0;
    sorted.forEach((id, k) => {
      if (k > 0 && cmp(keyOf(sorted[k - 1]), keyOf(id)) !== 0) rank = k;
      r.set(id, rank);
    });
    return r;
  };
  let ranks = toRanks((id) => inv.get(id));
  let classes = new Set(ranks.values()).size;
  for (let it = 0; it < ids.length; it++) {
    const cur = ranks;
    const next = toRanks((id) => [cur.get(id), ...mol.bondsOf(id)
      .map((b) => cur.get(mol.otherAtom(b, id)) * 4 + b.order).sort((a, b) => a - b)]);
    const n = new Set(next.values()).size;
    ranks = next;
    if (n === classes) break;
    classes = n;
  }
  // Break ties deterministically by id.
  const sorted = [...ids].sort((x, y) => ranks.get(x) - ranks.get(y) || x - y);
  const final = new Map();
  sorted.forEach((id, k) => final.set(id, k));
  return final;
}

function cmp(a, b) {
  for (let i = 0; i < Math.max(a.length, b.length); i++) {
    const x = a[i] ?? -Infinity, y = b[i] ?? -Infinity;
    if (x !== y) return x < y ? -1 : 1;
  }
  return 0;
}

function writeFragment(mol, start, ranks, aroAtoms, aroBonds, visited) {
  // Pass 1: DFS to find tree edges and ring closures.
  const children = new Map();
  const closuresAt = new Map(); // atom -> [{bond, partner}]
  const seenBond = new Set();
  const order = [];
  const sortedNbs = (id) => mol.bondsOf(id).slice()
    .sort((x, y) => ranks.get(mol.otherAtom(x, id)) - ranks.get(mol.otherAtom(y, id)));
  const dfs = (root) => {
    const stack = [[root, null]];
    // Iterative DFS mirroring recursion order.
    const iters = new Map();
    visited.add(root); order.push(root); children.set(root, []);
    while (stack.length) {
      const [v] = stack[stack.length - 1];
      if (!iters.has(v)) iters.set(v, { list: sortedNbs(v), k: 0 });
      const it = iters.get(v);
      if (it.k >= it.list.length) { stack.pop(); continue; }
      const b = it.list[it.k++];
      if (seenBond.has(b.id)) continue;
      seenBond.add(b.id);
      const w = mol.otherAtom(b, v);
      if (visited.has(w)) {
        // Ring closure: opened at ancestor w, closed at v.
        if (!closuresAt.has(w)) closuresAt.set(w, []);
        if (!closuresAt.has(v)) closuresAt.set(v, []);
        const rc = { bond: b, open: w, close: v, digit: null };
        closuresAt.get(w).push(rc);
        closuresAt.get(v).push(rc);
      } else {
        visited.add(w); order.push(w);
        children.get(v).push({ atom: w, bond: b });
        children.set(w, []);
        stack.push([w, b]);
      }
    }
  };
  dfs(start);

  // Pass 2: emit.
  const freeDigits = [];
  let nextDigit = 1;
  const takeDigit = () => {
    if (freeDigits.length) { freeDigits.sort((a, b) => a - b); return freeDigits.shift(); }
    return nextDigit++;
  };
  const digitStr = (d) => (d < 10 ? String(d) : `%${d}`);
  const bondStr = (b) => {
    if (aroBonds.has(b.id) && aroAtoms.has(b.a1) && aroAtoms.has(b.a2)) return '';
    if (b.order === 2) return '=';
    if (b.order === 3) return '#';
    if (b.order === 4) return '$';
    if (aroAtoms.has(b.a1) && aroAtoms.has(b.a2)) return '-';
    return '';
  };
  const out = [];
  const emit = (root) => {
    // Iterative emission using an explicit stack of tasks.
    const tasks = [{ type: 'atom', id: root }];
    while (tasks.length) {
      const t = tasks.pop();
      if (t.type === 'text') { out.push(t.text); continue; }
      const id = t.id;
      out.push(atomStr(mol, id, aroAtoms, aroBonds));
      const rcs = closuresAt.get(id) || [];
      // Closings first (digit already assigned), then openings.
      const closes = rcs.filter((rc) => rc.close === id && rc.digit !== null);
      const opens = rcs.filter((rc) => rc.open === id && rc.digit === null);
      for (const rc of closes) out.push(digitStr(rc.digit));
      for (const rc of opens) {
        rc.digit = takeDigit();
        out.push(bondStr(rc.bond) + digitStr(rc.digit));
      }
      // Release closed digits only now so they are not reused on this atom.
      for (const rc of closes) freeDigits.push(rc.digit);
      const ch = children.get(id) || [];
      // Push in reverse so that the first child is processed first.
      const seq = [];
      ch.forEach((c, k) => {
        const last = k === ch.length - 1;
        if (!last) seq.push({ type: 'text', text: '(' });
        seq.push({ type: 'text', text: bondStr(c.bond) });
        seq.push({ type: 'atom', id: c.atom });
        if (!last) seq.push({ type: 'branchEnd' });
      });
      // Branch ends must be emitted after the whole subtree; convert to text tasks.
      for (let k = seq.length - 1; k >= 0; k--) {
        const s = seq[k];
        tasks.push(s.type === 'branchEnd' ? { type: 'text', text: ')' } : s);
      }
    }
  };
  emit(start);
  return out.join('');
}

function atomStr(mol, id, aroAtoms, aroBonds) {
  const atom = mol.getAtom(id);
  const aro = aroAtoms.has(id);
  let organic;
  if (aro) {
    const bonds = mol.bondsOf(id);
    const isAro = (b) => aroBonds.has(b.id) && aroAtoms.has(b.a1) && aroAtoms.has(b.a2);
    const aroB = bonds.filter(isAro);
    let used = aroB.length;
    for (const b of bonds) if (!isAro(b)) used += b.order;
    organic = canWriteOrganic(mol, id, true, used, aroB.some((b) => b.order === 2));
  } else {
    organic = canWriteOrganic(mol, id, false);
  }
  const sym = aro ? atom.el.toLowerCase() : atom.el;
  if (organic) return sym;
  const h = mol.implicitH(id);
  let s = '[';
  if (atom.isotope) s += atom.isotope;
  s += sym;
  if (h) s += h === 1 ? 'H' : `H${h}`;
  if (atom.charge) {
    const sign = atom.charge > 0 ? '+' : '-';
    s += Math.abs(atom.charge) === 1 ? sign : `${sign}${Math.abs(atom.charge)}`;
  }
  if (atom.mapNo) s += `:${atom.mapNo}`;
  return `${s}]`;
}
