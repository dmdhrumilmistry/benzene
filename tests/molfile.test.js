import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readMolfile, writeMolfile } from '../js/core/molfile.js';
import { parseSmiles } from '../js/core/smiles.js';
import { BOND_STEREO } from '../js/core/molecule.js';

function formula(mol) {
  const c = {};
  for (const a of mol.atomList()) {
    c[a.el] = (c[a.el] || 0) + 1;
    const h = mol.implicitH(a.id);
    if (h) c.H = (c.H || 0) + h;
  }
  return Object.keys(c).sort().map((k) => k + c[k]).join('');
}

const ETHANOL = `ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
`;

test('reads V2000, flips y and normalizes bond length', () => {
  const m = readMolfile(ETHANOL);
  assert.equal(formula(m), 'C2H6O1');
  const [a, b] = m.atomList();
  assert.ok(Math.abs(Math.hypot(a.x - b.x, a.y - b.y) - 1) < 1e-3);
  assert.ok(b.y < a.y, 'y flipped (atom 2 was higher in y-up space)');
});

test('reads SDF first record, charges, isotopes, radicals', () => {
  const text = `acetate
  x

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   2  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
M  CHG  1   3  -1
M  ISO  1   1  13
M  END
> <NAME>
acetate

$$$$
second
`;
  const m = readMolfile(text);
  const atoms = m.atomList();
  assert.equal(atoms[2].charge, -1);
  assert.equal(atoms[0].isotope, 13);
  assert.equal(atoms[1].isotope, 14);
  assert.equal(formula(m), 'C2H3O2');
  const rad = readMolfile(ETHANOL.replace('M  END', 'M  RAD  1   1   2'));
  assert.equal(rad.atomList()[0].radical, 1);
  assert.equal(formula(rad), 'C2H5O1');
});

test('reads aromatic bond type 4 and kekulizes', () => {
  const smi = parseSmiles('c1ccccc1');
  let text = writeMolfile(smi);
  text = text.replace(/^(\s+\d+\s+\d+)\s+[12](\s+0\s+0\s+0\s+0)$/gm, (s, a, b) => `${a}  4${b}`);
  assert.equal((text.match(/  4  0  0  0  0$/gm) || []).length, 6);
  const m = readMolfile(text);
  assert.equal(m.bondList().filter((b) => b.order === 2).length, 3);
  assert.equal(formula(m), 'C6H6');
});

test('round-trips through writeMolfile', () => {
  for (const smi of ['CC(=O)Oc1ccccc1C(=O)O', '[NH4+]', '[13CH4]', 'C[CH2]', '[O-]C(=O)C', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C']) {
    const m = parseSmiles(smi);
    const text = writeMolfile(m, { name: smi });
    assert.ok(text.endsWith('M  END'));
    assert.equal(text.split('\n')[0], smi);
    const r = readMolfile(text);
    assert.equal(formula(r), formula(m), smi);
    assert.equal(r.bondCount, m.bondCount);
  }
});

test('writes stereo codes and scaled, flipped coordinates', () => {
  const m = parseSmiles('CC(O)N');
  const b = m.bondList()[1];
  b.stereo = BOND_STEREO.WEDGE;
  m.bondList()[2].stereo = BOND_STEREO.HASH;
  const text = writeMolfile(m);
  const lines = text.split('\n');
  assert.match(lines[3], /^\s+4\s+3 .*V2000$/);
  const bondLines = lines.slice(8, 11);
  assert.equal(bondLines[1].slice(9, 12).trim(), '1');
  assert.equal(bondLines[2].slice(9, 12).trim(), '6');
  const r = readMolfile(text);
  assert.equal(r.bondList()[1].stereo, BOND_STEREO.WEDGE);
  assert.equal(r.bondList()[2].stereo, BOND_STEREO.HASH);
  const x1 = parseFloat(lines[4].slice(0, 10)), y1 = parseFloat(lines[4].slice(10, 20));
  const x2 = parseFloat(lines[5].slice(0, 10)), y2 = parseFloat(lines[5].slice(10, 20));
  assert.ok(Math.abs(Math.hypot(x1 - x2, y1 - y2) - 1.5) < 1e-3);
  const a1 = m.atomList()[0], a2 = m.atomList()[1];
  assert.ok(Math.sign(y2 - y1) === -Math.sign(a2.y - a1.y) || Math.abs(a2.y - a1.y) < 1e-6);
});

test('lays out molecules without coordinates', () => {
  const text = ETHANOL.replace(/^\s+[\d.]+\s+[\d.]+\s+0\.0000/gm, '    0.0000    0.0000    0.0000');
  const m = readMolfile(text);
  const [a, b] = m.atomList();
  assert.ok(Math.abs(Math.hypot(a.x - b.x, a.y - b.y) - 1) < 1e-3);
});

test('reads basic V3000', () => {
  const text = `
  test

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1.5 0 0 0
M  V30 3 N 2.25 1.3 0 0 CHG=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 2 2 3
M  V30 END BOND
M  V30 END CTAB
M  END`;
  const m = readMolfile(text);
  assert.equal(m.atomCount, 3);
  assert.equal(m.atomList()[2].charge, 1);
  assert.equal(m.bondList()[0].stereo, BOND_STEREO.WEDGE);
  assert.equal(formula(m), 'C2H6N1');
});

test('throws on garbage', () => {
  assert.throws(() => readMolfile('hello'), /Invalid molfile/);
});
