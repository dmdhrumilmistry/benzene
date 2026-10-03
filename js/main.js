// Benzene application bootstrap: wires the editor to menus, tool palette,
// dialogs, side panels, storage and keyboard shortcuts.

import { Editor } from './editor/editor.js';
import { ChemDocument, normalizeMolecule, stripExplicitHydrogens } from './editor/document.js';
import { exportSVG, DEFAULT_STYLE } from './editor/renderer.js';
import { STRUCTURE_LIBRARY, localNameLookup } from './editor/templates.js';
import { parseSmiles } from './core/smiles.js';
import { readMolfile, writeMolfile } from './core/molfile.js';
import { layoutMolecule } from './core/layout.js';
import { computeProperties } from './core/properties.js';
import { ELEMENT_LIST } from './core/elements.js';
import { TOOL_GROUPS, ICONS, iconFor, matchToolId } from './ui/tools.js';
import { PropertiesPanel, NmrPanel, LibraryPanel, smilesFor, thumbnailSVG } from './ui/panels.js';
import { Library, Autosave, Settings } from './ui/storage.js';
import { download, safeFilename, exportPNG, parseStructureText, readFileAsDocument } from './ui/fileio.js';
import { nameToMolfile } from './ui/pubchem.js';

const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
const isMac = /Mac|iPhone|iPad/.test(navigator.platform);
const MOD = isMac ? '⌘' : 'Ctrl+';

const ACS_STYLE = { bondLength: 19.2, lineWidth: 0.8, fontSize: 13.3, doubleGap: 0.18 };
const SETTINGS_DEFAULTS = { showAtomNumbers: false, showCarbons: false, colorAtoms: true, grid: false, libraryEntryId: null };

class App {
  constructor() {
    this.settings = Settings.load(SETTINGS_DEFAULTS);
    this.editor = new Editor($('#canvas'), { parseSmiles, layoutMolecule });
    this.editor.setStyle({
      showAtomNumbers: this.settings.showAtomNumbers, showCarbons: this.settings.showCarbons, colorAtoms: this.settings.colorAtoms,
    });
    $('#canvas-wrap').classList.toggle('grid', this.settings.grid);
    this.currentEntryId = null;
    this.props = new PropertiesPanel($('#tab-props'), this);
    this.nmr = new NmrPanel($('#tab-nmr'), this);
    this.library = new LibraryPanel($('#tab-library'), this);

    this._buildMenus();
    this._buildToolbox();
    this._buildElementBar();
    this._buildQuickActions();
    this._buildPeriodicTable();
    this._bindDialogs();
    this._bindTabs();
    this._bindKeyboard();
    this._bindDragDrop();
    this._bindEditorEvents();
    this._restore();
    this.editor.setTool({ name: 'bond', order: 1, stereo: 'none' });
    this._refresh();
  }

  // ================================================================ helpers
  toast(msg, { error = false, ms = 2200 } = {}) {
    const t = $('#toast');
    t.textContent = msg;
    t.classList.toggle('error', error);
    t.classList.add('show');
    clearTimeout(this._toastT);
    this._toastT = setTimeout(() => t.classList.remove('show'), ms);
  }

  docName() { return $('#doc-name').value.trim() || 'Untitled'; }
  setDocName(n) { $('#doc-name').value = n; this._autosave(); }

  status(msg) { $('#status-msg').textContent = msg; }

  _refresh() {
    const ed = this.editor;
    $('#empty-hint').classList.toggle('hidden', !ed.doc.isEmpty());
    this.props.update();
    this.nmr.update();
    try {
      const p = ed.doc.mol.isEmpty() ? null : computeProperties(ed.doc.mol);
      $('#status-formula').innerHTML = p ? `${p.formulaHTML} · ${p.molecularWeight.toFixed(2)}` : '';
    } catch { $('#status-formula').textContent = ''; }
    $('#qa-undo').disabled = !ed.undoStack.length;
    $('#qa-redo').disabled = !ed.redoStack.length;
  }

  _autosave() {
    clearTimeout(this._saveT);
    this._saveT = setTimeout(() => Autosave.save(this.editor.doc.toJSON(), this.docName()), 400);
  }

  _restore() {
    const params = new URLSearchParams(location.search);
    const smi = params.get('smiles');
    if (smi) {
      try { this.loadDocFromMol(parseSmiles(smi), params.get('name') || 'Untitled'); return; } catch (e) { this.toast(`Bad SMILES in URL: ${e.message}`, { error: true }); }
    }
    const saved = Autosave.load();
    if (saved?.doc) {
      try {
        this.editor.setDocument(ChemDocument.fromJSON(saved.doc));
        $('#doc-name').value = saved.name || 'Untitled';
      } catch (e) { console.warn('Autosave restore failed', e); }
    }
  }

  loadDocFromMol(mol, name = 'Untitled') {
    const doc = new ChemDocument();
    doc.mol = mol;
    normalizeMolecule(doc.mol);
    this.editor.setDocument(doc);
    this.setDocName(name);
    this.currentEntryId = null;
  }

  confirmDiscard() {
    if (this.editor.doc.isEmpty() || !this.editor.undoStack.length) return true;
    return confirm('Discard the current drawing? (Save it to the library first if you want to keep it.)');
  }

  // ================================================================== menus
  _menuDefs() {
    const ed = this.editor;
    const s = this.settings;
    return [
      {
        label: 'File', items: [
          { label: 'New', key: `${MOD}N`, action: () => this.newDocument() },
          { label: 'Open file…', key: `${MOD}O`, action: () => $('#file-input').click() },
          { label: 'Save to library', key: `${MOD}S`, action: () => this.saveToLibrary() },
          { label: 'Download document (.benzene.json)', action: () => this.downloadNative() },
          '-',
          { label: 'Import SMILES / MOL…', action: () => this.openSmilesDialog() },
          { label: 'Name to structure…', action: () => this.openNameDialog() },
          '-',
          { label: 'Export image (PNG / SVG)…', key: `${MOD}E`, action: () => this.openExportDialog() },
          { label: 'Export MOL file', action: () => this.exportMol() },
          { label: 'Export SMILES file', action: () => this.exportSmiles() },
          '-',
          { label: 'Copy SMILES', action: () => this.copyText(smilesFor(ed.doc.mol, this._selOrNull()), 'SMILES copied') },
          { label: 'Copy MOL block', action: () => this.copyText(writeMolfile(this._selMol(), { name: this.docName() }), 'MOL block copied') },
          { label: 'Copy image', action: () => this.copyImage() },
          { label: 'Copy shareable link', action: () => this.copyLink() },
        ],
      },
      {
        label: 'Edit', items: [
          { label: 'Undo', key: `${MOD}Z`, action: () => ed.undo() },
          { label: 'Redo', key: isMac ? '⇧⌘Z' : 'Ctrl+Y', action: () => ed.redo() },
          '-',
          { label: 'Cut', key: `${MOD}X`, action: () => this.cut() },
          { label: 'Copy', key: `${MOD}C`, action: () => this.copy() },
          { label: 'Paste', key: `${MOD}V`, action: () => this.paste() },
          { label: 'Duplicate', key: `${MOD}D`, action: () => { ed.copySelection(); ed.paste(); } },
          { label: 'Delete', key: 'Del', action: () => ed.deleteSelection() },
          '-',
          { label: 'Select all', key: `${MOD}A`, action: () => ed.selectAll() },
          { label: 'Deselect', key: 'Esc', action: () => ed.clearSelection() },
        ],
      },
      {
        label: 'View', items: [
          { label: 'Zoom in', key: `${MOD}+`, action: () => ed.zoomBy(1.25) },
          { label: 'Zoom out', key: `${MOD}−`, action: () => ed.zoomBy(0.8) },
          { label: 'Fit to window', key: `${MOD}0`, action: () => ed.fitToView() },
          '-',
          { label: 'Show atom numbers', check: () => s.showAtomNumbers, action: () => this.toggleSetting('showAtomNumbers') },
          { label: 'Show carbon labels', check: () => s.showCarbons, action: () => this.toggleSetting('showCarbons') },
          { label: 'Color atom labels', check: () => s.colorAtoms, action: () => this.toggleSetting('colorAtoms') },
          { label: 'Show grid', check: () => s.grid, action: () => { s.grid = !s.grid; Settings.save(s); $('#canvas-wrap').classList.toggle('grid', s.grid); } },
          '-',
          { label: 'Toggle dark mode', action: () => this.toggleTheme() },
        ],
      },
      {
        label: 'Structure', items: [
          { label: 'Clean up structure', key: `${MOD}⇧K`, action: () => ed.cleanStructure() },
          { label: 'Rotate 90° clockwise', action: () => ed.rotateSelection(Math.PI / 2) },
          { label: 'Rotate 90° counter-clockwise', action: () => ed.rotateSelection(-Math.PI / 2) },
          { label: 'Rotate 15°', action: () => ed.rotateSelection(Math.PI / 12) },
          { label: 'Flip horizontal', action: () => ed.flipSelection('h') },
          { label: 'Flip vertical', action: () => ed.flipSelection('v') },
          '-',
          { label: 'Add explicit hydrogens', action: () => this.addExplicitH() },
          { label: 'Remove explicit hydrogens', action: () => ed.change(() => stripExplicitHydrogens(ed.doc.mol)) },
          { label: 'Check structure', action: () => this.checkStructure() },
          '-',
          { label: 'Templates…', key: 'Shift+T', action: () => this.openTemplates() },
          { label: 'Periodic table…', action: () => $('#dlg-periodic').showModal() },
        ],
      },
      {
        label: 'Help', items: [
          { label: 'Keyboard shortcuts', key: '?', action: () => this.openHelp() },
          { label: 'About Benzene', action: () => $('#dlg-about').showModal() },
          { label: 'Report an issue', action: () => window.open('https://github.com/dmdhrumilmistry/benzene/issues', '_blank', 'noopener') },
        ],
      },
    ];
  }

  _buildMenus() {
    const bar = $('#menubar');
    bar.innerHTML = '';
    const defs = this._menuDefs();
    for (const def of defs) {
      const menu = document.createElement('div');
      menu.className = 'menu';
      const btn = document.createElement('button');
      btn.textContent = def.label;
      const list = document.createElement('div');
      list.className = 'menu-list';
      menu.append(btn, list);
      const fill = () => {
        list.innerHTML = '';
        for (const item of def.items) {
          if (item === '-') { const sep = document.createElement('div'); sep.className = 'menu-sep'; list.appendChild(sep); continue; }
          const it = document.createElement('button');
          it.className = 'menu-item';
          const check = item.check ? (item.check() ? '✓' : '') : '';
          it.innerHTML = `<span class="check">${check}</span><span>${item.label}</span>${item.key ? `<span class="shortcut">${item.key}</span>` : ''}`;
          it.addEventListener('click', () => { this._closeMenus(); item.action(); });
          list.appendChild(it);
        }
      };
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const open = menu.classList.contains('open');
        this._closeMenus();
        if (!open) { fill(); menu.classList.add('open'); }
      });
      btn.addEventListener('mouseenter', () => {
        if ($$('.menu.open').length && !menu.classList.contains('open')) { this._closeMenus(); fill(); menu.classList.add('open'); }
      });
      bar.appendChild(menu);
    }
    document.addEventListener('click', () => this._closeMenus());
  }

  _closeMenus() { $$('.menu.open').forEach((m) => m.classList.remove('open')); }

  // ================================================================ toolbox
  _buildToolbox() {
    const box = $('#toolbox');
    box.innerHTML = '';
    for (const g of TOOL_GROUPS) {
      const group = document.createElement('div');
      group.className = 'tool-group';
      group.innerHTML = `<div class="tool-group-title">${g.title}</div>`;
      for (const t of g.tools) {
        const b = document.createElement('button');
        b.className = 'tool';
        b.dataset.tool = t.id;
        b.title = t.title;
        b.setAttribute('aria-label', t.title);
        b.innerHTML = iconFor(t);
        b.addEventListener('click', () => this.selectTool(t));
        group.appendChild(b);
      }
      box.appendChild(group);
    }
    const g = document.createElement('div');
    g.className = 'tool-group';
    g.innerHTML = '<div class="tool-group-title">Library</div>';
    const tb = document.createElement('button');
    tb.className = 'tool';
    tb.title = 'Structure templates (Shift+T)';
    tb.innerHTML = ICONS.templates;
    tb.addEventListener('click', () => this.openTemplates());
    g.appendChild(tb);
    box.appendChild(g);
    this._updateAtomToolLabel();
  }

  selectTool(def) {
    this.editor.setTool(def.tool);
  }

  _syncToolUI() {
    const id = matchToolId(this.editor.tool);
    $$('.tool[data-tool]').forEach((b) => b.classList.toggle('active', b.dataset.tool === id));
    $$('.el-btn[data-el]').forEach((b) => b.classList.toggle('active', this.editor.tool.name === 'atom' && b.dataset.el === this.editor.element));
    const t = this.editor.tool;
    const hints = {
      select: 'Click or drag to select · drag selection to move · drag the round handle to rotate · double-click selects a fragment',
      erase: 'Click or drag over atoms, bonds and objects to delete them',
      bond: 'Click empty space or an atom to add a bond · drag to set direction · click a bond to change it',
      chain: 'Drag to draw a zigzag chain',
      ring: 'Click empty space, an atom (spiro) or a bond (fused) to add the ring',
      atom: `Click atoms to change them to ${this.editor.element} · click empty space to place ${this.editor.element}`,
      charge: 'Click an atom to change its charge',
      text: 'Click empty space to add text (use _2 for subscript, ^+ for superscript) · click an atom to type a label like OMe, CO2H, NH2',
      arrow: 'Drag to draw an arrow (hold Alt for free angle)',
      pan: 'Drag to pan the canvas',
    };
    this.status(hints[t.name] || 'Ready');
  }

  _buildElementBar() {
    const bar = $('#elementbar');
    const els = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si'];
    bar.innerHTML = els.map((e) => `<button class="el-btn" data-el="${e}" title="Atom tool: ${e}">${e}</button>`).join('')
      + '<button class="el-btn" id="btn-periodic" title="Periodic table">⋯</button>';
    $$('.el-btn[data-el]', bar).forEach((b) => b.addEventListener('click', () => this.selectElement(b.dataset.el)));
    $('#btn-periodic').addEventListener('click', () => $('#dlg-periodic').showModal());
  }

  selectElement(el) {
    this.editor.element = el;
    // If atoms are selected, change them directly.
    const ids = [...this.editor.selection.atoms];
    if (ids.length && this.editor.tool.name === 'select') {
      this.editor.change(() => ids.forEach((id) => this.editor.setAtomElement(id, el)));
      return;
    }
    this._updateAtomToolLabel();
    this.editor.setTool({ name: 'atom' });
  }

  _updateAtomToolLabel() {
    const l = document.getElementById('atom-tool-label');
    if (l) {
      l.textContent = this.editor.element;
      l.setAttribute('font-size', this.editor.element.length > 1 ? '12' : '15');
    }
  }

  _buildQuickActions() {
    const ed = this.editor;
    const acts = [
      ['qa-undo', 'undo', `Undo (${MOD}Z)`, () => ed.undo()],
      ['qa-redo', 'redo', 'Redo', () => ed.redo()],
      '-',
      ['qa-clean', 'clean', 'Clean up structure', () => ed.cleanStructure()],
      ['qa-rotate', 'rotate', 'Rotate 15° (selection or all)', () => ed.rotateSelection(Math.PI / 12)],
      ['qa-fliph', 'flipH', 'Flip horizontal', () => ed.flipSelection('h')],
      ['qa-flipv', 'flipV', 'Flip vertical', () => ed.flipSelection('v')],
      ['qa-delete', 'trash', 'Delete selection', () => ed.deleteSelection()],
      '-',
      ['qa-zoomout', 'zoomOut', 'Zoom out', () => ed.zoomBy(0.8)],
      ['qa-zoomin', 'zoomIn', 'Zoom in', () => ed.zoomBy(1.25)],
      ['qa-fit', 'fit', 'Fit to window', () => ed.fitToView()],
      '-',
      ['qa-export', 'export', 'Export image', () => this.openExportDialog()],
    ];
    const host = $('#quick-actions');
    for (const a of acts) {
      if (a === '-') { const s = document.createElement('span'); s.className = 'act-sep'; host.appendChild(s); continue; }
      const [id, icon, title, fn] = a;
      const b = document.createElement('button');
      b.className = 'act-btn';
      b.id = id;
      b.title = title;
      b.setAttribute('aria-label', title);
      b.innerHTML = ICONS[icon];
      b.addEventListener('click', fn);
      host.appendChild(b);
    }
  }

  _buildPeriodicTable() {
    const host = $('#periodic-table');
    const cells = [];
    for (const el of ELEMENT_LIST) {
      const row = el.row >= 9 ? el.row + 1 : el.row; // leave a gap row before f-block
      cells.push(`<button class="pt-cell" style="grid-row:${row};grid-column:${el.col};color:${el.color}" data-el="${el.symbol}" title="${el.name} (${el.z}) — ${el.mass}"><b>${el.symbol}</b><span>${el.z}</span></button>`);
    }
    host.innerHTML = cells.join('') + '<div style="grid-row:9;grid-column:1/-1;height:6px"></div>';
    $$('.pt-cell', host).forEach((b) => b.addEventListener('click', () => {
      $('#dlg-periodic').close();
      this.selectElement(b.dataset.el);
    }));
  }

  // ================================================================ dialogs
  _bindDialogs() {
    $$('dialog [data-close]').forEach((b) => b.addEventListener('click', () => b.closest('dialog').close()));
    $$('dialog').forEach((d) => d.addEventListener('keydown', (e) => e.stopPropagation()));
    $$('[data-action="import-smiles"]').forEach((b) => b.addEventListener('click', () => this.openSmilesDialog()));
    $$('[data-action="name-lookup"]').forEach((b) => b.addEventListener('click', () => this.openNameDialog()));

    $('#smiles-ok').addEventListener('click', () => this._importSmiles());
    $('#smiles-input').addEventListener('keydown', (e) => { if (e.key === 'Enter' && !e.shiftKey && !$('#smiles-input').value.includes('\n')) { e.preventDefault(); this._importSmiles(); } });

    $('#name-ok').addEventListener('click', () => this._importName());
    $('#name-input').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); this._importName(); } });

    for (const id of ['export-format', 'export-scale', 'export-bg', 'export-style']) $(`#${id}`).addEventListener('change', () => this._updateExportPreview());
    $('#export-ok').addEventListener('click', () => this._doExport());
    $('#export-copy').addEventListener('click', () => this.copyImage());

    $('#template-search').addEventListener('input', () => this._renderTemplates());

    $('#file-input').addEventListener('change', async (e) => {
      const f = e.target.files[0];
      e.target.value = '';
      if (f) this.openFile(f);
    });
    $('#library-input').addEventListener('change', async (e) => {
      const f = e.target.files[0];
      e.target.value = '';
      if (!f) return;
      try { const n = Library.importJSON(await f.text()); this.toast(`Imported ${n} structure${n === 1 ? '' : 's'}`); this.library.render(); } catch (err) { this.toast(err.message, { error: true }); }
    });

    $('#btn-theme').addEventListener('click', () => this.toggleTheme());
    $('#doc-name').addEventListener('change', () => this._autosave());
    $('#doc-name').addEventListener('keydown', (e) => { e.stopPropagation(); if (e.key === 'Enter') e.target.blur(); });
  }

  openSmilesDialog() {
    $('#smiles-error').textContent = '';
    $('#dlg-smiles').showModal();
    $('#smiles-input').focus();
    $('#smiles-input').select();
  }

  _importSmiles() {
    const text = $('#smiles-input').value;
    try {
      const doc = parseStructureText(text);
      this.editor.insertMolecule(doc.mol);
      $('#dlg-smiles').close();
      this.editor.fitToView();
    } catch (e) {
      $('#smiles-error').textContent = e.message;
    }
  }

  openNameDialog() {
    $('#name-error').textContent = '';
    $('#dlg-name').showModal();
    $('#name-input').focus();
    $('#name-input').select();
  }

  async _importName() {
    const name = $('#name-input').value.trim();
    if (!name) return;
    const err = $('#name-error');
    err.textContent = '';
    try {
      let mol;
      const local = localNameLookup(name);
      if (local) mol = parseSmiles(local);
      else {
        err.textContent = 'Searching PubChem…';
        err.style.color = 'var(--muted)';
        const molText = await nameToMolfile(name);
        mol = readMolfile(molText);
        stripExplicitHydrogens(mol);
      }
      err.style.color = '';
      const wasEmpty = this.editor.doc.isEmpty();
      this.editor.insertMolecule(mol);
      if (wasEmpty) this.setDocName(name.charAt(0).toUpperCase() + name.slice(1));
      this.editor.fitToView();
      $('#dlg-name').close();
    } catch (e) {
      err.style.color = '';
      err.textContent = e.message.includes('fetch') ? 'Could not reach PubChem (offline?)' : e.message;
    }
  }

  openTemplates() {
    $('#template-search').value = '';
    this._renderTemplates();
    $('#dlg-templates').showModal();
  }

  _renderTemplates() {
    const q = $('#template-search').value.trim().toLowerCase();
    const host = $('#template-grid');
    this._thumbCache = this._thumbCache || new Map();
    let html = '';
    for (const group of STRUCTURE_LIBRARY) {
      const items = group.items.filter(([n]) => !q || n.toLowerCase().includes(q));
      if (!items.length) continue;
      html += `<h5>${group.category}</h5><div class="template-cat">`;
      for (const [name, smi] of items) {
        let thumb = this._thumbCache.get(smi);
        if (thumb === undefined) {
          try {
            const d = new ChemDocument();
            d.mol = parseSmiles(smi);
            normalizeMolecule(d.mol);
            thumb = thumbnailSVG(d, this.editor.style);
          } catch (e) { thumb = `<span class="muted small">${e.message}</span>`; }
          this._thumbCache.set(smi, thumb);
        }
        html += `<button class="template-card" data-smiles="${smi}" data-name="${name}"><div class="thumb">${thumb}</div><span>${name}</span></button>`;
      }
      html += '</div>';
    }
    host.innerHTML = html || '<div class="placeholder">No matches</div>';
    $$('.template-card', host).forEach((b) => b.addEventListener('click', () => {
      try {
        const wasEmpty = this.editor.doc.isEmpty();
        this.editor.insertMolecule(parseSmiles(b.dataset.smiles));
        if (wasEmpty && this.docName() === 'Untitled') this.setDocName(b.dataset.name);
        $('#dlg-templates').close();
        this.editor.fitToView();
      } catch (e) { this.toast(e.message, { error: true }); }
    }));
  }

  openHelp() {
    const sections = [
      ['Tools', [['Select', 'V'], ['Eraser', 'E'], ['Single / double / triple bond', '1 / 2 / 3'], ['Wedge / hash bond', 'W / Q'], ['Chain', 'A'], ['Benzene ring', 'R'], ['Ring of size n', '3 … 8'], ['Text', 'T'], ['Reaction arrow', 'X'], ['Templates', 'Shift+T'], ['Pan', 'Space + drag / right-drag']]],
      ['Hover an atom', [['Change element', 'C N O S P F H I'], ['Chlorine / Bromine', 'L / B'], ['Increase / decrease charge', '+ / −'], ['Add a bond', '1 / 2 / 3'], ['Type a label (OMe, CO2H …)', 'Enter'], ['Delete', 'Del']]],
      ['Hover a bond', [['Set order', '1 / 2 / 3'], ['Wedge / hash (again = flip)', 'W / H'], ['Wavy', 'Y'], ['Plain', '0'], ['Delete', 'Del']]],
      ['Edit', [['Undo / redo', `${MOD}Z / ${isMac ? '⇧⌘Z' : 'Ctrl+Y'}`], ['Cut / copy / paste', `${MOD}X / C / V`], ['Duplicate', `${MOD}D`], ['Select all', `${MOD}A`], ['Nudge selection', 'Arrow keys (Shift = more)'], ['Clean up structure', `${MOD}Shift+K`]]],
      ['File & view', [['New', `${MOD}N`], ['Open', `${MOD}O`], ['Save to library', `${MOD}S`], ['Export image', `${MOD}E`], ['Zoom', `${MOD}+ / ${MOD}− / wheel+${isMac ? '⌘' : 'Ctrl'}`], ['Fit to window', `${MOD}0`]]],
      ['Tips', [['Paste SMILES or MOL text directly', MOD + 'V'], ['Drop .mol/.sdf/.smi/.svg/.png files', 'Drag & drop'], ['Exported PNG/SVG stay editable', 'Re-open them'], ['Share a structure', 'File → Copy shareable link']]],
    ];
    $('#help-content').innerHTML = sections.map(([title, rows]) => `<h5>${title}</h5>${rows.map(([a, k]) => `<div class="help-row"><span>${a}</span><kbd>${k}</kbd></div>`).join('')}`).join('');
    $('#dlg-help').showModal();
  }

  // ================================================================= export
  _exportStyle() {
    const base = { ...this.editor.style, background: null };
    return $('#export-style').value === 'acs' ? { ...base, ...ACS_STYLE } : base;
  }

  /** Document to export: selection if any, else everything. */
  _exportDoc() {
    const ed = this.editor;
    if (!ed.hasSelection()) return ed.doc;
    const d = new ChemDocument();
    d.mol = ed.doc.mol.subMolecule(ed.selectedAtomIds());
    for (const o of ed.doc.objects) if (ed.selection.objects.has(o.id)) d.addObject({ ...o });
    return d;
  }

  openExportDialog() {
    if (this.editor.doc.isEmpty()) { this.toast('Nothing to export yet'); return; }
    this._updateExportPreview();
    $('#dlg-export').showModal();
  }

  _updateExportPreview() {
    const bg = $('#export-bg').value;
    const { svg } = exportSVG(this._exportDoc(), this._exportStyle(), { embed: false, background: bg === 'transparent' ? null : bg });
    $('#export-preview').innerHTML = svg.replace(/^<\?xml[^>]*>/, '');
  }

  async _doExport() {
    const fmt = $('#export-format').value;
    const bg = $('#export-bg').value;
    const embed = $('#export-embed').checked;
    const doc = this._exportDoc();
    const style = this._exportStyle();
    try {
      if (fmt === 'svg') {
        const { svg } = exportSVG(doc, style, { embed, background: bg === 'transparent' ? null : bg, scale: Number($('#export-scale').value) });
        download(safeFilename(this.docName(), 'svg'), svg, 'image/svg+xml');
      } else {
        const blob = await exportPNG(doc, style, { scale: Number($('#export-scale').value), background: bg, embed });
        download(safeFilename(this.docName(), 'png'), blob);
      }
      $('#dlg-export').close();
    } catch (e) { this.toast(`Export failed: ${e.message}`, { error: true }); }
  }

  async copyImage() {
    if (this.editor.doc.isEmpty()) return;
    try {
      const blob = await exportPNG(this._exportDoc(), this._exportStyle(), { scale: 2, background: '#ffffff', embed: true });
      await navigator.clipboard.write([new ClipboardItem({ 'image/png': blob })]);
      this.toast('Image copied to clipboard');
    } catch (e) { this.toast(`Copy failed: ${e.message}`, { error: true }); }
  }

  _selOrNull() {
    const ids = this.editor.selectedAtomIds();
    return ids.size ? [...ids] : null;
  }

  _selMol() {
    const ids = this._selOrNull();
    return ids ? this.editor.doc.mol.subMolecule(ids) : this.editor.doc.mol;
  }

  exportMol() {
    if (this.editor.doc.mol.isEmpty()) { this.toast('Nothing to export yet'); return; }
    download(safeFilename(this.docName(), 'mol'), writeMolfile(this._selMol(), { name: this.docName() }), 'chemical/x-mdl-molfile');
  }

  exportSmiles() {
    const smi = smilesFor(this.editor.doc.mol, this._selOrNull());
    if (!smi) { this.toast('Nothing to export yet'); return; }
    download(safeFilename(this.docName(), 'smi'), `${smi} ${this.docName()}\n`, 'chemical/x-daylight-smiles');
  }

  downloadNative() {
    download(safeFilename(this.docName(), 'benzene.json'), JSON.stringify({ ...this.editor.doc.toJSON(), name: this.docName() }, null, 1), 'application/json');
  }

  copyText(text, msg) {
    if (!text) { this.toast('Nothing to copy'); return; }
    navigator.clipboard?.writeText(text).then(() => this.toast(msg), () => this.toast('Clipboard unavailable', { error: true }));
  }

  copyLink() {
    const smi = smilesFor(this.editor.doc.mol);
    if (!smi) { this.toast('Nothing to share'); return; }
    const url = new URL(location.href);
    url.search = '';
    url.searchParams.set('smiles', smi);
    if (this.docName() !== 'Untitled') url.searchParams.set('name', this.docName());
    this.copyText(url.toString(), 'Link copied');
  }

  // ================================================================== files
  newDocument() {
    if (!this.confirmDiscard()) return;
    this.editor.setDocument(new ChemDocument());
    this.setDocName('Untitled');
    this.currentEntryId = null;
  }

  async openFile(file) {
    try {
      const { doc, name } = await readFileAsDocument(file);
      if (this.editor.doc.isEmpty()) {
        this.editor.setDocument(doc);
        this.setDocName(name);
        this.currentEntryId = null;
      } else {
        // Add to the current drawing.
        this.editor.insertMolecule(doc.mol);
        if (doc.objects.length) this.editor.change(() => doc.objects.forEach((o) => this.editor.doc.addObject({ ...o, id: undefined })));
        this.editor.fitToView();
      }
      this.toast(`Opened ${file.name}`);
    } catch (e) {
      this.toast(`Could not open ${file.name}: ${e.message}`, { error: true, ms: 4000 });
    }
  }

  saveToLibrary(name = null) {
    const ed = this.editor;
    if (ed.doc.isEmpty()) { this.toast('Nothing to save yet'); return; }
    const n = (name ?? this.docName()).trim() || 'Untitled';
    this.setDocName(n);
    let formula = '';
    try { formula = computeProperties(ed.doc.mol)?.formula || ''; } catch { /* ignore */ }
    const payload = { name: n, doc: ed.doc.toJSON(), thumbnail: thumbnailSVG(ed.doc, ed.style), smiles: smilesFor(ed.doc.mol) || '', formula };
    try {
      const existing = this.currentEntryId ? Library.get(this.currentEntryId) : null;
      if (existing && existing.name === n) Library.update(existing.id, payload);
      else this.currentEntryId = Library.add(payload).id;
      this.toast(`Saved “${n}” to library`);
      this.library.render();
    } catch (e) { this.toast(e.message, { error: true, ms: 4000 }); }
  }

  openLibraryEntry(id) {
    const e = Library.get(id);
    if (!e || !this.confirmDiscard()) return;
    this.editor.setDocument(ChemDocument.fromJSON(e.doc));
    this.setDocName(e.name);
    this.currentEntryId = id;
  }

  insertLibraryEntry(id) {
    const e = Library.get(id);
    if (!e) return;
    const doc = ChemDocument.fromJSON(e.doc);
    if (!doc.mol.isEmpty()) this.editor.insertMolecule(doc.mol);
    this.editor.fitToView();
  }

  // ============================================================ clipboard
  copy() {
    const r = this.editor.copySelection();
    if (!r) return;
    let text = '';
    try { text = r.mol.isEmpty() ? '' : smilesFor(r.mol); } catch { /* ignore */ }
    this._lastCopiedText = text;
    if (text) navigator.clipboard?.writeText(text).catch(() => {});
    this.toast('Copied');
  }

  cut() { this.copy(); this.editor.deleteSelection(); }

  async paste() {
    // Prefer system clipboard text when it differs from what we copied (e.g. SMILES from elsewhere).
    let text = '';
    try { text = (await navigator.clipboard?.readText()) || ''; } catch { /* permission denied */ }
    if (text.trim() && text.trim() !== this._lastCopiedText) {
      try {
        const doc = parseStructureText(text);
        this.editor.insertMolecule(doc.mol, this.editor.pointer);
        return;
      } catch { /* not a structure: fall back to internal clipboard */ }
    }
    if (this.editor.clipboard) this.editor.paste();
  }

  // ================================================================ misc
  toggleSetting(key) {
    this.settings[key] = !this.settings[key];
    Settings.save(this.settings);
    this.editor.setStyle({ [key]: this.settings[key] });
  }

  toggleTheme() {
    const root = document.documentElement;
    const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
    try { localStorage.setItem('benzene.theme', root.dataset.theme); } catch { /* ignore */ }
  }

  addExplicitH() {
    const ed = this.editor;
    const ids = ed.selectedAtomIds();
    ed.change(() => {
      const mol = ed.doc.mol;
      for (const a of mol.atomList()) {
        if (ids.size && !ids.has(a.id)) continue;
        if (a.el === 'H') continue;
        const n = mol.implicitH(a.id);
        for (let i = 0; i < n; i++) {
          const h = ed.addBondFromAtom(a.id, 1, 'none', 'H');
          h.el = 'H';
          if (a.hCount !== null) a.hCount = Math.max(0, a.hCount - 1);
        }
      }
    });
  }

  checkStructure() {
    const mol = this.editor.doc.mol;
    const bad = [...mol.atoms.keys()].filter((id) => mol.hasValenceError(id));
    if (!bad.length) { this.toast('No valence problems found ✓'); return; }
    this.editor.selection = { atoms: new Set(bad), bonds: new Set(), objects: new Set() };
    this.editor.setTool({ name: 'select' });
    this.toast(`${bad.length} atom${bad.length > 1 ? 's have' : ' has'} an unusual valence (selected)`, { error: true });
  }

  // ================================================================ events
  _bindTabs() {
    $$('.tab').forEach((t) => t.addEventListener('click', () => {
      $$('.tab').forEach((x) => x.classList.toggle('active', x === t));
      $$('.tab-panel').forEach((p) => p.classList.toggle('active', p.id === `tab-${t.dataset.tab}`));
      this.nmr.setActive(t.dataset.tab === 'nmr');
      if (t.dataset.tab === 'library') this.library.render();
    }));
  }

  _bindEditorEvents() {
    const ed = this.editor;
    ed.addEventListener('change', () => { this._refresh(); this._autosave(); });
    ed.addEventListener('selection', () => this.props.update());
    ed.addEventListener('tool', () => this._syncToolUI());
    ed.addEventListener('view', () => { $('#status-zoom').textContent = `${Math.round(ed.view.zoom * 100)}%`; });
    ed.addEventListener('pointer', (e) => { $('#status-pos').textContent = `${e.detail.x.toFixed(1)}, ${e.detail.y.toFixed(1)}`; });
    ed.addEventListener('hover', (e) => {
      const h = e.detail;
      if (!h) { this._syncToolUI(); return; }
      const mol = ed.doc.mol;
      if (h.type === 'atom') {
        const a = mol.getAtom(h.id);
        const hs = mol.implicitH(h.id);
        this.status(`Atom ${a.label || a.el}${hs ? `H${hs > 1 ? hs : ''}` : ''}${a.charge ? ` (charge ${a.charge > 0 ? '+' : ''}${a.charge})` : ''} — type an element key, +/−, or Enter for a label`);
      } else if (h.type === 'bond') {
        const b = mol.getBond(h.id);
        this.status(`${['', 'Single', 'Double', 'Triple'][b.order] || ''} bond${b.stereo !== 'none' ? ` (${b.stereo})` : ''} — press 1/2/3, W, H, Y or Del`);
      }
    });
    $('#status-zoom').textContent = '100%';
  }

  _bindKeyboard() {
    const ed = this.editor;
    window.addEventListener('keydown', (e) => {
      const tag = document.activeElement?.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || document.querySelector('dialog[open]')) return;
      const mod = e.ctrlKey || e.metaKey;
      const k = e.key.toLowerCase();
      if (mod) {
        const map = {
          z: () => (e.shiftKey ? ed.redo() : ed.undo()), y: () => ed.redo(), c: () => this.copy(), x: () => this.cut(),
          v: () => this.paste(), a: () => ed.selectAll(), s: () => this.saveToLibrary(), o: () => $('#file-input').click(),
          e: () => this.openExportDialog(), n: () => this.newDocument(), d: () => { ed.copySelection(); ed.paste(); },
          '0': () => ed.fitToView(), '=': () => ed.zoomBy(1.25), '+': () => ed.zoomBy(1.25), '-': () => ed.zoomBy(0.8),
          k: () => (e.shiftKey ? ed.cleanStructure() : null),
        };
        if (map[k]) { e.preventDefault(); map[k](); }
        return;
      }
      if (e.altKey) return;
      if (ed.handleKey(e)) { e.preventDefault(); return; }
      if (e.key === '?') { this.openHelp(); return; }
      if (e.key === 'T' && e.shiftKey) { this.openTemplates(); return; }
      // Tool shortcuts (only when nothing is hovered).
      for (const g of TOOL_GROUPS) {
        for (const t of g.tools) {
          if (t.key && t.key === e.key) { e.preventDefault(); this.selectTool(t); return; }
        }
      }
    });
    window.addEventListener('keyup', (e) => ed.handleKeyUp(e));
  }

  _bindDragDrop() {
    const wrap = $('#canvas-wrap');
    wrap.addEventListener('dragover', (e) => { e.preventDefault(); wrap.classList.add('dragover'); });
    wrap.addEventListener('dragleave', () => wrap.classList.remove('dragover'));
    wrap.addEventListener('drop', (e) => {
      e.preventDefault();
      wrap.classList.remove('dragover');
      const f = e.dataTransfer.files[0];
      if (f) this.openFile(f);
      else {
        const text = e.dataTransfer.getData('text/plain');
        if (text) {
          try { this.editor.insertMolecule(parseStructureText(text).mol, this.editor.toModel(e.clientX, e.clientY)); } catch (err) { this.toast(err.message, { error: true }); }
        }
      }
    });
    window.addEventListener('beforeunload', () => Autosave.save(this.editor.doc.toJSON(), this.docName()));
  }
}

export { DEFAULT_STYLE };

window.benzene = new App();
