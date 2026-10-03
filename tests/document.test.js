import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ChemDocument, stripExplicitHydrogens, normalizeMolecule } from '../js/editor/document.js';
import { renderDocument, exportSVG, extractDocFromSVG, atomLabelInfo, DEFAULT_STYLE } from '../js/editor/renderer.js';
import { parseSmiles } from '../js/core/smiles.js';
import { Molecule } from '../js/core/molecule.js';

test('document JSON round trip keeps molecule and objects', () => {
  const doc = new ChemDocument();
  doc.mol = parseSmiles('CCO');
  doc.addObject({ type: 'arrow', style: 'forward', x1: 0, y1: 0, x2: 2, y2: 0 });
  doc.addObject({ type: 'text', x: 1, y: 1, text: 'H_2O', size: 1 });
  const copy = ChemDocument.fromJSON(JSON.parse(JSON.stringify(doc.toJSON())));
  assert.equal(copy.mol.atomCount, 3);
  assert.equal(copy.mol.bondCount, 2);
  assert.equal(copy.objects.length, 2);
  assert.equal(copy.objects[1].text, 'H_2O');
});

test('SVG export embeds an editable document', () => {
  const doc = new ChemDocument();
  doc.mol = parseSmiles('c1ccccc1O');
  const { svg, width, height } = exportSVG(doc);
  assert.ok(svg.startsWith('<?xml'));
  assert.ok(width > 0 && height > 0);
  const json = extractDocFromSVG(svg);
  assert.equal(ChemDocument.fromJSON(json).mol.atomCount, 7);
});

test('renderer labels heteroatoms with implicit hydrogens', () => {
  const doc = new ChemDocument();
  doc.mol = parseSmiles('CCO');
  const markup = renderDocument(doc, DEFAULT_STYLE);
  assert.match(markup, />O</);
  assert.match(markup, />H</);
  const o = doc.mol.atomList().find((a) => a.el === 'O');
  const info = atomLabelInfo(doc.mol, o, DEFAULT_STYLE);
  assert.equal(info.hCount, 1);
  // Plain carbons are not labelled.
  const c = doc.mol.atomList().find((a) => a.el === 'C');
  assert.equal(atomLabelInfo(doc.mol, c, DEFAULT_STYLE), null);
});

test('stripExplicitHydrogens converts H atoms to implicit', () => {
  const m = new Molecule();
  const c = m.addAtom({ el: 'C' });
  for (let i = 0; i < 4; i++) m.addBond(c.id, m.addAtom({ el: 'H', x: i }).id);
  stripExplicitHydrogens(m);
  assert.equal(m.atomCount, 1);
  assert.equal(m.implicitH(c.id), 4);
});

test('normalizeMolecule scales to unit bond length and centres', () => {
  const m = new Molecule();
  const a = m.addAtom({ x: 10, y: 10 });
  const b = m.addAtom({ x: 13, y: 10 });
  m.addBond(a.id, b.id);
  normalizeMolecule(m, 0, 0);
  assert.ok(Math.abs(m.averageBondLength() - 1) < 1e-9);
  const box = m.bbox();
  assert.ok(Math.abs(box.cx) < 1e-9 && Math.abs(box.cy) < 1e-9);
});
