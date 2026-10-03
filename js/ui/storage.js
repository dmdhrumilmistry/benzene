// Browser storage: the structure library, autosave and user settings.
// Everything lives in localStorage on the user's machine.

const KEYS = {
  library: 'benzene.library.v1',
  autosave: 'benzene.autosave.v1',
  settings: 'benzene.settings.v1',
};

function read(key, fallback) {
  try {
    const v = localStorage.getItem(key);
    return v ? JSON.parse(v) : fallback;
  } catch {
    return fallback;
  }
}

function write(key, value) {
  try {
    localStorage.setItem(key, JSON.stringify(value));
    return true;
  } catch (e) {
    console.warn('Storage write failed', e);
    return false;
  }
}

const uid = () => Date.now().toString(36) + Math.random().toString(36).slice(2, 8);

export const Library = {
  list() {
    return read(KEYS.library, []);
  },
  /** Save a new entry; returns the entry or throws if storage is full. */
  add({ name, doc, thumbnail, smiles = '', formula = '' }) {
    const list = this.list();
    const now = new Date().toISOString();
    const entry = { id: uid(), name: name || 'Untitled', doc, thumbnail, smiles, formula, created: now, updated: now };
    list.unshift(entry);
    if (!write(KEYS.library, list)) throw new Error('Browser storage is full — export your library and delete some entries.');
    return entry;
  },
  update(id, patch) {
    const list = this.list();
    const i = list.findIndex((e) => e.id === id);
    if (i < 0) return null;
    list[i] = { ...list[i], ...patch, updated: new Date().toISOString() };
    if (!write(KEYS.library, list)) throw new Error('Browser storage is full.');
    return list[i];
  },
  remove(id) {
    write(KEYS.library, this.list().filter((e) => e.id !== id));
  },
  get(id) {
    return this.list().find((e) => e.id === id) || null;
  },
  findByName(name) {
    return this.list().find((e) => e.name === name) || null;
  },
  exportJSON() {
    return JSON.stringify({ format: 'benzene-library', version: 1, entries: this.list() }, null, 2);
  },
  /** Merge entries from an exported library file. Returns number imported. */
  importJSON(text) {
    const json = JSON.parse(text);
    const entries = json.entries || (Array.isArray(json) ? json : null);
    if (!entries) throw new Error('Not a Benzene library file');
    const list = this.list();
    const ids = new Set(list.map((e) => e.id));
    let n = 0;
    for (const e of entries) {
      if (!e.doc) continue;
      if (ids.has(e.id)) e.id = uid();
      list.push(e);
      n++;
    }
    if (!write(KEYS.library, list)) throw new Error('Browser storage is full.');
    return n;
  },
};

export const Autosave = {
  save(doc, name) { write(KEYS.autosave, { doc, name, time: Date.now() }); },
  load() { return read(KEYS.autosave, null); },
  clear() { try { localStorage.removeItem(KEYS.autosave); } catch { /* ignore */ } },
};

export const Settings = {
  load(defaults) { return { ...defaults, ...read(KEYS.settings, {}) }; },
  save(settings) { write(KEYS.settings, settings); },
};
