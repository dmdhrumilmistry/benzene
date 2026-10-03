// A Benzene document: one Molecule (possibly several fragments) plus
// graphical objects (reaction arrows, curved arrows, text, plus signs).

import { Molecule } from '../core/molecule.js';

export const FORMAT_VERSION = 1;

export class ChemDocument {
  constructor() {
    this.mol = new Molecule();
    /** @type {Array<object>} arrows: {id,type:'arrow',style,x1,y1,x2,y2,curve}; text: {id,type:'text',x,y,text,size} */
    this.objects = [];
    this.nextObjectId = 1;
  }

  addObject(obj) {
    const o = { ...obj, id: obj.id ?? this.nextObjectId++ };
    if (o.id >= this.nextObjectId) this.nextObjectId = o.id + 1;
    this.objects.push(o);
    return o;
  }

  getObject(id) { return this.objects.find((o) => o.id === id); }

  removeObject(id) { this.objects = this.objects.filter((o) => o.id !== id); }

  isEmpty() { return this.mol.isEmpty() && this.objects.length === 0; }

  /** Bounding box in model units of everything in the document. */
  bbox() {
    let b = this.mol.bbox();
    const grow = (x, y) => {
      if (!b) b = { minX: x, minY: y, maxX: x, maxY: y };
      b.minX = Math.min(b.minX, x); b.minY = Math.min(b.minY, y);
      b.maxX = Math.max(b.maxX, x); b.maxY = Math.max(b.maxY, y);
    };
    for (const o of this.objects) {
      if (o.type === 'arrow') { grow(o.x1, o.y1); grow(o.x2, o.y2); }
      else if (o.type === 'text') {
        const lines = String(o.text).split('\n');
        const w = Math.max(...lines.map((l) => l.length)) * 0.28 * (o.size || 1);
        grow(o.x, o.y - 0.4); grow(o.x + w, o.y + lines.length * 0.45 * (o.size || 1));
      }
    }
    if (!b) return null;
    b = { ...b };
    b.width = b.maxX - b.minX; b.height = b.maxY - b.minY;
    b.cx = (b.minX + b.maxX) / 2; b.cy = (b.minY + b.maxY) / 2;
    return b;
  }

  clone() {
    return ChemDocument.fromJSON(this.toJSON());
  }

  toJSON() {
    return {
      format: 'benzene',
      version: FORMAT_VERSION,
      molecule: this.mol.toJSON(),
      objects: this.objects.map((o) => ({ ...o })),
    };
  }

  static fromJSON(json) {
    const doc = new ChemDocument();
    if (!json) return doc;
    doc.mol = Molecule.fromJSON(json.molecule || { atoms: [], bonds: [] });
    for (const o of json.objects || []) doc.addObject(o);
    return doc;
  }
}

/**
 * Remove explicit hydrogen atoms attached to a single heavy atom by a plain
 * single bond (they become implicit). Used when importing structures from
 * sources such as PubChem that list every H explicitly.
 */
export function stripExplicitHydrogens(mol) {
  for (const atom of mol.atomList()) {
    if (atom.el !== 'H' || atom.isotope || atom.charge) continue;
    const bonds = mol.bondsOf(atom.id);
    if (bonds.length !== 1 || bonds[0].order !== 1 || bonds[0].stereo !== 'none') continue;
    const heavy = mol.getAtom(mol.otherAtom(bonds[0], atom.id));
    if (heavy.el === 'H') continue;
    // Keep the heavy atom's H count unchanged if it had an explicit count.
    if (heavy.hCount !== null) heavy.hCount += 1;
    mol.removeAtom(atom.id);
  }
  return mol;
}

/** Normalise imported molecules to bond length 1 and centre them on (cx, cy). */
export function normalizeMolecule(mol, cx = 0, cy = 0) {
  if (mol.isEmpty()) return mol;
  const l = mol.averageBondLength();
  if (mol.bondCount && Math.abs(l - 1) > 0.02) mol.scale(1 / l);
  const b = mol.bbox();
  mol.translate(cx - b.cx, cy - b.cy);
  return mol;
}
