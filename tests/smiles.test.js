import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseSmiles, writeSmiles, kekulize } from '../js/core/smiles.js';
import { Molecule } from '../js/core/molecule.js';

function formula(mol) {
  const c = {};
  for (const a of mol.atomList()) {
    c[a.el] = (c[a.el] || 0) + 1;
    const h = mol.implicitH(a.id);
    if (h) c.H = (c.H || 0) + h;
  }
  const keys = Object.keys(c).sort((a, b) => (a === 'C' ? -1 : b === 'C' ? 1 : a === 'H' ? -1 : b === 'H' ? 1 : a.localeCompare(b)));
  return keys.map((k) => k + (c[k] > 1 ? c[k] : '')).join('');
}

/** Connectivity fingerprint: refined atom classes + bond order multiset. */
function graphKey(mol) {
  let cls = new Map(mol.atomList().map((a) => [a.id, `${a.el}${a.charge}h${mol.implicitH(a.id)}d${mol.degree(a.id)}`]));
  for (let i = 0; i < 6; i++) {
    const next = new Map();
    for (const a of mol.atomList()) {
      const nb = mol.bondsOf(a.id).map((b) => `${b.order}${cls.get(mol.otherAtom(b, a.id))}`).sort().join(',');
      next.set(a.id, `${cls.get(a.id)}[${nb}]`);
    }
    cls = next;
  }
  const bonds = mol.bondList().map((b) => b.order).sort().join('');
  return `${[...cls.values()].sort().join('|')}#${bonds}`;
}

const SET = [
  ['CCO', 'C2H6O'], ['c1ccccc1', 'C6H6'], ['CC(=O)Oc1ccccc1C(=O)O', 'C9H8O4'],
  ['Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C8H10N4O2'], ['c1ccc2ccccc2c1', 'C10H8'], ['C1CC1', 'C3H6'],
  ['[NH4+]', 'H4N'], ['[O-]C(=O)C', 'C2H3O2'], ['OC[C@H](N)C(=O)O', 'C3H7NO3'], ['c1cc[nH]c1', 'C4H5N'],
  ['c1ccoc1', 'C4H4O'], ['c1ccsc1', 'C4H4S'], ['n1ccccc1', 'C5H5N'], ['C#N', 'CHN'], ['[13CH4]', 'CH4'],
  ['CC(C)(C)C', 'C5H12'], ['C1CCC2CCCCC2C1', 'C10H18'], ['O=C1NC(=O)c2ccccc21', 'C8H5NO2'],
  ['c1ccc2[nH]ccc2c1', 'C8H7N'], ['Clc1ccc(Br)cc1', 'C6H4BrCl'], ['[Na+].[Cl-]', 'ClNa'],
  ['C=C/C=C/C', 'C5H8'], ['N#Cc1ccccc1', 'C7H5N'], ['CS(=O)(=O)C', 'C2H6O2S'],
];

test('parses formulas correctly', () => {
  for (const [smi, f] of SET) assert.equal(formula(parseSmiles(smi)), f, smi);
});

test('round-trips formula and connectivity (aromatic and Kekulé output)', () => {
  for (const [smi] of SET) {
    const m1 = parseSmiles(smi);
    for (const aromatic of [true, false]) {
      const out = writeSmiles(m1, { aromatic });
      const m2 = parseSmiles(out);
      assert.equal(formula(m2), formula(m1), `${smi} -> ${out}`);
      assert.equal(graphKey(m2), graphKey(m1), `${smi} -> ${out}`);
      assert.equal(writeSmiles(m1, { aromatic }), out, 'deterministic');
    }
  }
});

test('kekulizes aromatic rings to alternating bonds', () => {
  const m = parseSmiles('c1ccccc1');
  const orders = m.bondList().map((b) => b.order).sort();
  assert.deepEqual(orders, [1, 1, 1, 2, 2, 2]);
  for (const a of m.atomList()) {
    assert.equal(a.hCount, null);
    assert.equal(m.implicitH(a.id), 1);
    assert.equal(m.bondsOf(a.id).filter((b) => b.order === 2).length, 1);
  }
  const ind = parseSmiles('c1ccc2[nH]ccc2c1');
  assert.equal(ind.bondList().filter((b) => b.order === 2).length, 4);
  const pyrene = parseSmiles('c1cc2ccc3cccc4ccc(c1)c2c34');
  assert.equal(pyrene.bondList().filter((b) => b.order === 2).length, 8);
});

test('writes lowercase aromatic SMILES and brackets where needed', () => {
  assert.equal(writeSmiles(parseSmiles('C1=CC=CC=C1')), 'c1ccccc1');
  assert.match(writeSmiles(parseSmiles('c1cc[nH]c1')), /\[nH\]/);
  assert.equal(writeSmiles(parseSmiles('[NH4+]')), '[NH4+]');
  assert.equal(writeSmiles(parseSmiles('[13CH4]')), '[13CH4]');
  assert.equal(writeSmiles(parseSmiles('[Na+].[Cl-]')), '[Na+].[Cl-]');
  assert.equal(writeSmiles(parseSmiles('[CH3]')), '[CH3]');
  assert.equal(writeSmiles(new Molecule()), '');
});

test('parses bracket atom details', () => {
  const m = parseSmiles('[13CH3:7][C@@H](F)[O-]');
  const [c13, cs, , o] = m.atomList();
  assert.equal(c13.isotope, 13);
  assert.equal(c13.hCount, 3);
  assert.equal(c13.mapNo, 7);
  assert.equal(cs.chirality, '@@');
  assert.equal(o.charge, -1);
  assert.equal(parseSmiles('[Fe+2]').atomList()[0].charge, 2);
  assert.equal(parseSmiles('[Fe++]').atomList()[0].charge, 2);
  assert.equal(parseSmiles('[se]1cccc1').atomList()[0].el, 'Se');
});

test('ring closures with %nn and bond symbols', () => {
  const m = parseSmiles('C%10CCCCC%10');
  assert.equal(m.bondCount, 6);
  const m2 = parseSmiles('C=1CCCCC1');
  assert.equal(m2.bondList().filter((b) => b.order === 2).length, 1);
  assert.equal(formula(parseSmiles('CCO ethanol')), 'C2H6O');
});

test('generates 2D coordinates', () => {
  const m = parseSmiles('CCO');
  const [a, b] = m.atomList();
  assert.ok(Math.abs(Math.hypot(a.x - b.x, a.y - b.y) - 1) < 1e-6);
});

test('throws helpful errors on invalid input', () => {
  for (const bad of ['C(C', 'C)C', 'C1CC', 'Xx', '[Zz]', 'C==C', 'c1cccc1', '[C', 'C%1', '=C']) {
    assert.throws(() => parseSmiles(bad), /Invalid SMILES/, bad);
  }
});

test('kekulize helper is exported and tolerant', () => {
  const m = new Molecule();
  const ids = [];
  for (let i = 0; i < 6; i++) ids.push(m.addAtom({ el: 'C' }).id);
  const bonds = ids.map((id, i) => m.addBond(id, ids[(i + 1) % 6], 1).id);
  assert.equal(kekulize(m, bonds).success, true);
  assert.equal(m.bondList().filter((b) => b.order === 2).length, 3);
});
