// File import/export helpers: downloads, SVG -> PNG rasterisation and
// embedding the editable document inside PNG files (tEXt chunk).

import { ChemDocument, normalizeMolecule, stripExplicitHydrogens } from '../editor/document.js';
import { exportSVG, extractDocFromSVG } from '../editor/renderer.js';
import { parseSmiles } from '../core/smiles.js';
import { readMolfile } from '../core/molfile.js';

const PNG_KEYWORD = 'benzene-document';

export function download(filename, data, mime = 'application/octet-stream') {
  const blob = data instanceof Blob ? data : new Blob([data], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

export function safeFilename(name, ext) {
  const base = (name || 'structure').trim().replace(/[^\w\-. ]+/g, '_').replace(/\s+/g, '_').slice(0, 80) || 'structure';
  return `${base}.${ext}`;
}

/** Rasterise an SVG string to a PNG Blob. */
export function svgToPngBlob(svgText, width, height, { scale = 2, background = null } = {}) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' }));
    img.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = Math.max(1, Math.ceil(width * scale));
      canvas.height = Math.max(1, Math.ceil(height * scale));
      const ctx = canvas.getContext('2d');
      if (background && background !== 'transparent') { ctx.fillStyle = background; ctx.fillRect(0, 0, canvas.width, canvas.height); }
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Could not render SVG')); };
    img.src = url;
  });
}

// ------------------------------------------------------------ PNG metadata
const CRC_TABLE = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();

function crc32(bytes) {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

function utf8ToBase64(str) {
  const bytes = new TextEncoder().encode(str);
  let bin = '';
  for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(bin);
}

function base64ToUtf8(b64) {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new TextDecoder().decode(bytes);
}

/** Insert a tEXt chunk (keyword + base64 UTF-8 text) right after IHDR. */
export function embedPngText(buffer, text, keyword = PNG_KEYWORD) {
  const src = new Uint8Array(buffer);
  const payload = new TextEncoder().encode(`${keyword}\0${utf8ToBase64(text)}`);
  const chunk = new Uint8Array(12 + payload.length);
  const dv = new DataView(chunk.buffer);
  dv.setUint32(0, payload.length);
  chunk.set([0x74, 0x45, 0x58, 0x74], 4); // 'tEXt'
  chunk.set(payload, 8);
  dv.setUint32(8 + payload.length, crc32(chunk.subarray(4, 8 + payload.length)));
  const ihdrEnd = 8 + 25; // signature + IHDR chunk (13 data bytes + 12)
  const out = new Uint8Array(src.length + chunk.length);
  out.set(src.subarray(0, ihdrEnd), 0);
  out.set(chunk, ihdrEnd);
  out.set(src.subarray(ihdrEnd), ihdrEnd + chunk.length);
  return out;
}

/** Read our tEXt chunk back from a PNG (or null). */
export function extractPngText(buffer, keyword = PNG_KEYWORD) {
  const bytes = new Uint8Array(buffer);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let off = 8;
  while (off + 12 <= bytes.length) {
    const length = dv.getUint32(off);
    const type = String.fromCharCode(...bytes.subarray(off + 4, off + 8));
    if (type === 'tEXt') {
      const data = bytes.subarray(off + 8, off + 8 + length);
      const nul = data.indexOf(0);
      const key = new TextDecoder('latin1').decode(data.subarray(0, nul));
      if (key === keyword) return base64ToUtf8(new TextDecoder('latin1').decode(data.subarray(nul + 1)));
    }
    if (type === 'IEND') break;
    off += 12 + length;
  }
  return null;
}

/** Export the document as PNG Blob (optionally embedding the editable document). */
export async function exportPNG(doc, style, { scale = 2, background = '#ffffff', embed = true } = {}) {
  const { svg, width, height } = exportSVG(doc, style, { embed: false, background: null });
  const blob = await svgToPngBlob(svg, width, height, { scale, background });
  if (!embed) return blob;
  const buf = await blob.arrayBuffer();
  return new Blob([embedPngText(buf, JSON.stringify(doc.toJSON()))], { type: 'image/png' });
}

// ------------------------------------------------------------------ import
/** Detect and parse structure text (SMILES, MOL/SDF, Benzene JSON, SVG). Returns a ChemDocument. */
export function parseStructureText(text, hint = '') {
  const t = text.trim();
  if (!t) throw new Error('Nothing to import');
  if (t.startsWith('{')) {
    const json = JSON.parse(t);
    if (json.format === 'benzene' || json.molecule) return ChemDocument.fromJSON(json);
    throw new Error('Unrecognised JSON file');
  }
  if (t.startsWith('<') && t.includes('<svg')) {
    const json = extractDocFromSVG(t);
    if (!json) throw new Error('This SVG has no embedded Benzene structure');
    return ChemDocument.fromJSON(json);
  }
  const doc = new ChemDocument();
  if (/M {2}END|V2000|V3000/.test(t) || /\.(mol|sdf|sd)$/i.test(hint)) {
    doc.mol = readMolfile(text);
    stripExplicitHydrogens(doc.mol);
  } else {
    // SMILES file: take the first token of the first non-empty line.
    const line = t.split(/\r?\n/).find((l) => l.trim());
    doc.mol = parseSmiles(line.trim().split(/\s+/)[0]);
  }
  normalizeMolecule(doc.mol);
  return doc;
}

/** Read a File into { doc, name }. */
export async function readFileAsDocument(file) {
  const name = file.name.replace(/\.[^.]+$/, '');
  if (/\.png$/i.test(file.name) || file.type === 'image/png') {
    const text = extractPngText(await file.arrayBuffer());
    if (!text) throw new Error('This PNG has no embedded Benzene structure');
    return { doc: ChemDocument.fromJSON(JSON.parse(text)), name };
  }
  const text = await file.text();
  return { doc: parseStructureText(text, file.name), name };
}
