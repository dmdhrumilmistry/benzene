import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseSmiles } from '../js/core/smiles.js';
import { layoutMolecule } from '../js/core/layout.js';

function metrics(mol) {
  const atoms = mol.atomList();
  let minNonBonded = Infinity, maxDev = 0;
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      if (mol.bondBetween(atoms[i].id, atoms[j].id)) continue;
      minNonBonded = Math.min(minNonBonded, Math.hypot(atoms[i].x - atoms[j].x, atoms[i].y - atoms[j].y));
    }
  }
  for (const b of mol.bondList()) {
    const p = mol.getAtom(b.a1), q = mol.getAtom(b.a2);
    maxDev = Math.max(maxDev, Math.abs(Math.hypot(p.x - q.x, p.y - q.y) - 1));
  }
  return { minNonBonded, maxDev };
}

const NON_BRIDGED = {
  aspirin: 'CC(=O)Oc1ccccc1C(=O)O',
  caffeine: 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
  naphthalene: 'c1ccc2ccccc2c1',
  steroid: 'CC(C)CCCC(C)C1CCC2C1(C)CCC1C2CC=C2CC(O)CCC12C',
  ibuprofen: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
  spiro: 'C1CCC2(CC1)CCCC2',
  pyrene: 'c1cc2ccc3cccc4ccc(c1)c2c34',
  hexamethylbenzene: 'Cc1c(C)c(C)c(C)c(C)c1C',
  ATP: 'Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O',
  decane: 'CCCCCCCCCC',
  peptide: 'NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O',
  macrocycle: 'C1CCCCCCCCCCC1',
};

for (const [name, smi] of Object.entries(NON_BRIDGED)) {
  test(`layout quality: ${name}`, () => {
    const m = parseSmiles(smi);
    const { minNonBonded, maxDev } = metrics(m);
    assert.ok(maxDev <= 0.05, `bond length deviation ${maxDev}`);
    assert.ok(minNonBonded > 0.5, `min non-bonded distance ${minNonBonded}`);
  });
}

test('benzene is a regular hexagon centred on the origin', () => {
  const m = parseSmiles('c1ccccc1');
  for (const a of m.atomList()) assert.ok(Math.abs(Math.hypot(a.x, a.y) - 1) < 1e-6);
});

test('triple bonds are linear', () => {
  const m = parseSmiles('CC#CC');
  const [a, b, c, d] = m.atomList();
  const ang = (p, q, r) => {
    const v1 = [p.x - q.x, p.y - q.y], v2 = [r.x - q.x, r.y - q.y];
    return Math.acos((v1[0] * v2[0] + v1[1] * v2[1]) / (Math.hypot(...v1) * Math.hypot(...v2)));
  };
  assert.ok(Math.abs(ang(a, b, c) - Math.PI) < 1e-6);
  assert.ok(Math.abs(ang(b, c, d) - Math.PI) < 1e-6);
});

test('chains zigzag at 120 degrees', () => {
  const m = parseSmiles('CCCC');
  const [a, , c] = m.atomList();
  assert.ok(Math.abs(Math.hypot(a.x - c.x, a.y - c.y) - Math.sqrt(3)) < 1e-6);
});

test('fragments laid out left to right without overlap', () => {
  const m = parseSmiles('c1ccccc1.[Na+].[Cl-]');
  const frags = m.fragments();
  const boxes = frags.map((f) => m.bbox(f));
  boxes.sort((p, q) => p.minX - q.minX);
  for (let i = 1; i < boxes.length; i++) assert.ok(boxes[i].minX > boxes[i - 1].maxX + 0.5);
  const bb = m.bbox();
  assert.ok(Math.abs(bb.cx) < 1e-6 && Math.abs(bb.cy) < 1e-6);
});

test('bridged systems get finite, non-coincident coordinates', () => {
  for (const smi of ['C1CC2CCC1C2', 'C1C2CC3CC1CC(C2)C3', 'C12C3C4C1C5C2C3C45']) {
    const m = parseSmiles(smi);
    layoutMolecule(m);
    for (const a of m.atomList()) assert.ok(Number.isFinite(a.x) && Number.isFinite(a.y));
    assert.ok(metrics(m).minNonBonded > 0.2, smi);
  }
});
