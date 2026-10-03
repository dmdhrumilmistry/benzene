import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseSmiles } from '../js/core/smiles.js';
import { computeProperties, isotopePattern, elementCounts } from '../js/core/properties.js';
import { analyzeRings } from '../js/core/rings.js';

const props = (smi) => computeProperties(parseSmiles(smi));
const near = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} expected ${b} ± ${tol}, got ${a}`);

test('aspirin formula, masses and descriptors', () => {
  const p = props('CC(=O)Oc1ccccc1C(=O)O');
  assert.equal(p.formula, 'C9H8O4');
  near(p.molecularWeight, 180.159, 0.01, 'MW');
  near(p.exactMass, 180.0423, 0.0005, 'exact mass');
  assert.equal(p.hBondDonors, 1);
  assert.equal(p.hBondAcceptors, 4);
  assert.equal(p.rings, 1);
  assert.equal(p.aromaticRings, 1);
  assert.equal(p.degreeOfUnsaturation, 6);
  near(p.tpsa, 63.6, 0.1, 'TPSA');
  assert.equal(p.rotatableBonds, 3);
  assert.equal(p.lipinskiViolations, 0);
});

test('elemental analysis sums to 100%', () => {
  const p = props('Cn1cnc2c1c(=O)n(C)c(=O)n2C');
  assert.equal(p.formula, 'C8H10N4O2');
  const total = p.elementalAnalysis.reduce((s, e) => s + e.percent, 0);
  near(total, 100, 1e-6);
  near(p.tpsa, 61.82, 0.1, 'caffeine TPSA');
});

test('charged species formula and m/z', () => {
  const p = props('C[N+](C)(C)C');
  assert.equal(p.formula, 'C4H12N+');
  near(p.mz, 74.0964, 0.001);
  assert.equal(props('CC(=O)[O-]').formula, 'C2H3O2-');
});

test('isotope pattern of chlorobenzene shows M+2 ≈ 32%', () => {
  const counts = elementCounts(parseSmiles('Clc1ccccc1'));
  const peaks = isotopePattern(counts);
  const m = peaks.find((p) => Math.abs(p.mass - 112.008) < 0.01);
  const m2 = peaks.find((p) => Math.abs(p.mass - 114.005) < 0.01);
  assert.ok(m && m2);
  near(m.abundance, 100, 0.01);
  near(m2.abundance, 32.0, 1.0, 'M+2');
});

test('ring perception: naphthalene and cubane', () => {
  const nap = analyzeRings(parseSmiles('c1ccc2ccccc2c1'));
  assert.equal(nap.rings.length, 2);
  assert.equal(nap.aromatic.atoms.size, 10);
  const cubane = analyzeRings(parseSmiles('C12C3C4C1C5C2C3C45'));
  assert.equal(cubane.rings.length, 5);
  assert.ok(cubane.rings.every((r) => r.length === 4));
  assert.equal(cubane.aromatic.atoms.size, 0);
});

test('heteroaromatic perception', () => {
  for (const smi of ['c1ccncc1', 'c1cc[nH]c1', 'c1ccoc1', 'c1ccsc1', 'c1ccc2[nH]ccc2c1']) {
    const m = parseSmiles(smi);
    const r = analyzeRings(m);
    assert.equal(r.aromatic.atoms.size, m.atomCount, smi);
  }
  assert.equal(analyzeRings(parseSmiles('C1=CCCCC1')).aromatic.atoms.size, 0);
});
