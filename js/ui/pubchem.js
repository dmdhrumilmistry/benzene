// Optional online lookups via the PubChem PUG REST API (CORS enabled).
// https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

const BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

/** Fetch a 2D SDF record for a compound name. Returns MOL text. */
export async function nameToMolfile(name) {
  const res = await fetch(`${BASE}/compound/name/${encodeURIComponent(name.trim())}/SDF?record_type=2d`);
  if (res.status === 404) throw new Error(`“${name}” was not found on PubChem`);
  if (!res.ok) throw new Error(`PubChem request failed (${res.status})`);
  const text = await res.text();
  return text.split('$$$$')[0];
}

/** Look up the IUPAC name and common title for a SMILES string. */
export async function smilesToName(smiles) {
  const body = new URLSearchParams({ smiles });
  const res = await fetch(`${BASE}/compound/smiles/property/IUPACName,Title/JSON`, { method: 'POST', body });
  if (res.status === 404) return null;
  if (!res.ok) throw new Error(`PubChem request failed (${res.status})`);
  const json = await res.json();
  const p = json?.PropertyTable?.Properties?.[0];
  if (!p || !p.CID) return null;
  return { cid: p.CID, iupac: p.IUPACName || null, title: p.Title || null };
}
