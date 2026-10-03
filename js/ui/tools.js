// Tool palette definitions and inline SVG icons.

import { RING_TEMPLATES } from '../editor/templates.js';

const svg = (inner) => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">${inner}</svg>`;

function polygon(n, r = 8.5, cx = 12, cy = 12.5) {
  const pts = [];
  const off = n === 4 || n === 8 ? Math.PI / n : 0;
  for (let i = 0; i < n; i++) {
    const t = -Math.PI / 2 + off + (i * 2 * Math.PI) / n;
    pts.push([cx + r * Math.cos(t), cy + r * Math.sin(t)]);
  }
  return pts;
}
function ringIcon(n, aromatic = false) {
  const pts = polygon(n);
  let inner = `<polygon points="${pts.map((p) => p.map((v) => v.toFixed(1)).join(',')).join(' ')}"/>`;
  if (aromatic) {
    const c = [12, 12.5];
    const k = 0.72;
    for (let i = n === 5 ? 1 : 0; i < n; i += 2) {
      if (n === 5 && i > 3) break;
      const a = pts[i], b = pts[(i + 1) % n];
      const ia = [c[0] + (a[0] - c[0]) * k, c[1] + (a[1] - c[1]) * k], ib = [c[0] + (b[0] - c[0]) * k, c[1] + (b[1] - c[1]) * k];
      inner += `<line x1="${ia[0].toFixed(1)}" y1="${ia[1].toFixed(1)}" x2="${ib[0].toFixed(1)}" y2="${ib[1].toFixed(1)}"/>`;
    }
  }
  return svg(inner);
}

export const ICONS = {
  select: svg('<rect x="3" y="3" width="13" height="13" stroke-dasharray="3 2.5"/><path d="M12 11l8 3-3.3 1.2L15.5 19z" fill="currentColor"/>'),
  erase: svg('<path d="M4 15l8.5-8.5a2 2 0 0 1 2.8 0l3.2 3.2a2 2 0 0 1 0 2.8L13 18H7z"/><path d="M9 10l5 5"/><path d="M13 20h7"/>'),
  pan: svg('<path d="M8 12V6a1.5 1.5 0 0 1 3 0v5M11 11V4.5a1.5 1.5 0 0 1 3 0V11M14 11V6a1.5 1.5 0 0 1 3 0v7c0 4-2.5 7-6 7-2.5 0-4-1.2-5.5-3.5L3.5 13a1.5 1.5 0 0 1 2.5-1.5L8 14"/>'),
  single: svg('<line x1="5" y1="19" x2="19" y2="5"/>'),
  double: svg('<line x1="4" y1="16.5" x2="16.5" y2="4"/><line x1="7.5" y1="20" x2="20" y2="7.5"/>'),
  triple: svg('<line x1="3" y1="15" x2="15" y2="3"/><line x1="5.5" y1="18.5" x2="18.5" y2="5.5"/><line x1="9" y1="21" x2="21" y2="9"/>'),
  wedge: svg('<polygon points="5,19 17.5,4.2 19.8,6.5" fill="currentColor"/>'),
  hash: svg('<g stroke-width="1.5"><line x1="6.5" y1="18.5" x2="5.5" y2="17.5"/><line x1="9.6" y1="16.2" x2="7.8" y2="14.4"/><line x1="12.7" y1="13.9" x2="10.1" y2="11.3"/><line x1="15.8" y1="11.6" x2="12.4" y2="8.2"/><line x1="18.9" y1="9.3" x2="14.7" y2="5.1"/></g>'),
  wavy: svg('<path d="M4 18c1-2 3-1 3.5-3s2.5-1 3.5-3 2.5-1 3.5-3 2.5-1 3.5-3"/>'),
  bold: svg('<line x1="5" y1="19" x2="19" y2="5" stroke-width="4.5" stroke-linecap="butt"/>'),
  dashed: svg('<line x1="5" y1="19" x2="19" y2="5" stroke-dasharray="2.5 2.5"/>'),
  chain: svg('<polyline points="3,16 8,8 13,16 18,8 22,14"/>'),
  atom: '<svg viewBox="0 0 24 24"><text x="12" y="17" text-anchor="middle" font-size="15" font-weight="700" fill="currentColor" font-family="Arial, sans-serif" id="atom-tool-label">N</text></svg>',
  chargePlus: svg('<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>'),
  chargeMinus: svg('<circle cx="12" cy="12" r="8"/><path d="M8 12h8"/>'),
  text: svg('<path d="M5 6V4h14v2M12 4v16M9 20h6"/>'),
  arrowForward: svg('<line x1="3" y1="12" x2="18" y2="12"/><path d="M15 8l5 4-5 4" fill="currentColor"/>'),
  arrowEquilibrium: svg('<path d="M3 9.5h16l-4-3.5"/><path d="M21 14.5H5l4 3.5"/>'),
  arrowResonance: svg('<line x1="6" y1="12" x2="18" y2="12"/><path d="M16 8.5l4 3.5-4 3.5zM8 8.5L4 12l4 3.5z" fill="currentColor"/>'),
  arrowRetro: svg('<path d="M3 10h14M3 14h14M14 6l6 6-6 6"/>'),
  arrowCurved: svg('<path d="M4 17C5 7 15 5 19 12"/><path d="M15.5 11.5l4 1 .5-4.2" fill="currentColor"/>'),
  arrowFishhook: svg('<path d="M4 17C5 7 15 5 19 12"/><path d="M19.5 12.5l.3-4.8-3.6 2.4z" fill="currentColor"/>'),
  arrowNoreaction: svg('<line x1="3" y1="12" x2="18" y2="12"/><path d="M15 8l5 4-5 4" fill="currentColor"/><path d="M8 8.5l4 7M12 8.5l-4 7"/>'),
  arrowDashed: svg('<line x1="3" y1="12" x2="17" y2="12" stroke-dasharray="3 2.5"/><path d="M15 8l5 4-5 4" fill="currentColor"/>'),
  templates: svg('<rect x="3" y="3" width="7.5" height="7.5" rx="1.5"/><rect x="13.5" y="3" width="7.5" height="7.5" rx="1.5"/><rect x="3" y="13.5" width="7.5" height="7.5" rx="1.5"/><path d="M17.2 14v6.5M14 17.2h6.5"/>'),
  undo: svg('<path d="M9 14L4 9l5-5"/><path d="M4 9h10.5a5.5 5.5 0 0 1 0 11H11"/>'),
  redo: svg('<path d="M15 14l5-5-5-5"/><path d="M20 9H9.5a5.5 5.5 0 0 0 0 11H13"/>'),
  clean: svg('<path d="M12 3l1.8 4.2L18 9l-4.2 1.8L12 15l-1.8-4.2L6 9l4.2-1.8z"/><path d="M18 15l.9 2.1L21 18l-2.1.9L18 21l-.9-2.1L15 18l2.1-.9z"/>'),
  rotate: svg('<path d="M20 12a8 8 0 1 1-2.3-5.7"/><path d="M20 4v5h-5"/>'),
  flipH: svg('<path d="M12 3v18" stroke-dasharray="2 2"/><path d="M9 7L3 17h6zM15 7l6 10h-6z"/>'),
  flipV: svg('<path d="M3 12h18" stroke-dasharray="2 2"/><path d="M7 9L17 3v6zM7 15l10 6v-6z"/>'),
  zoomIn: svg('<circle cx="10.5" cy="10.5" r="6.5"/><path d="M15.5 15.5L21 21M10.5 7.5v6M7.5 10.5h6"/>'),
  zoomOut: svg('<circle cx="10.5" cy="10.5" r="6.5"/><path d="M15.5 15.5L21 21M7.5 10.5h6"/>'),
  fit: svg('<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/>'),
  export: svg('<path d="M12 3v12M7 10l5 5 5-5"/><path d="M4 17v3h16v-3"/>'),
  trash: svg('<path d="M4 7h16M10 11v6M14 11v6M6 7l1 13h10l1-13M9 7V4h6v3"/>'),
};

export const TOOL_GROUPS = [
  {
    title: 'Edit',
    tools: [
      { id: 'select', icon: 'select', title: 'Select / move (V). Drag the handle to rotate', tool: { name: 'select' }, key: 'v' },
      { id: 'erase', icon: 'erase', title: 'Eraser (E)', tool: { name: 'erase' }, key: 'e' },
      { id: 'pan', icon: 'pan', title: 'Pan (hold Space)', tool: { name: 'pan' } },
      { id: 'text', icon: 'text', title: 'Text / atom label (T). Click an atom to type a label like OMe, CO2H, NH2', tool: { name: 'text' }, key: 't' },
    ],
  },
  {
    title: 'Bonds',
    tools: [
      { id: 'single', icon: 'single', title: 'Single bond (1). Click a bond to cycle order', tool: { name: 'bond', order: 1, stereo: 'none' }, key: '1' },
      { id: 'double', icon: 'double', title: 'Double bond (2). Click a double bond to move it', tool: { name: 'bond', order: 2, stereo: 'none' }, key: '2' },
      { id: 'triple', icon: 'triple', title: 'Triple bond (3)', tool: { name: 'bond', order: 3, stereo: 'none' }, key: '3' },
      { id: 'chain', icon: 'chain', title: 'Chain (A) - drag to draw a carbon chain', tool: { name: 'chain' }, key: 'a' },
      { id: 'wedge', icon: 'wedge', title: 'Wedge bond (W). Click again to flip', tool: { name: 'bond', order: 1, stereo: 'wedge' }, key: 'w' },
      { id: 'hash', icon: 'hash', title: 'Hashed wedge bond (Q). Click again to flip', tool: { name: 'bond', order: 1, stereo: 'hash' }, key: 'q' },
      { id: 'wavy', icon: 'wavy', title: 'Wavy bond (unknown stereo)', tool: { name: 'bond', order: 1, stereo: 'wavy' } },
      { id: 'bold', icon: 'bold', title: 'Bold bond', tool: { name: 'bond', order: 1, stereo: 'bold' } },
      { id: 'dashed', icon: 'dashed', title: 'Dashed bond (partial / H-bond)', tool: { name: 'bond', order: 1, stereo: 'dashed' } },
      { id: 'atom', icon: 'atom', title: 'Atom tool - click atoms to set the selected element', tool: { name: 'atom' } },
    ],
  },
  {
    title: 'Rings',
    tools: RING_TEMPLATES.map((t) => ({
      id: t.id, title: `${t.label}${t.id === 'benzene' ? ' (R)' : t.size >= 3 && !t.aromatic ? ` (${t.size})` : ''}`,
      iconHTML: ringIcon(t.size, t.aromatic), tool: { name: 'ring', template: t },
      key: t.id === 'benzene' ? 'r' : t.aromatic ? null : String(t.size),
    })),
  },
  {
    title: 'Charges',
    tools: [
      { id: 'plus', icon: 'chargePlus', title: 'Positive charge (click atom)', tool: { name: 'charge', delta: 1 }, key: '+' },
      { id: 'minus', icon: 'chargeMinus', title: 'Negative charge (click atom)', tool: { name: 'charge', delta: -1 }, key: '-' },
    ],
  },
  {
    title: 'Arrows',
    tools: [
      { id: 'arrow', icon: 'arrowForward', title: 'Reaction arrow (X)', tool: { name: 'arrow', style: 'forward' }, key: 'x' },
      { id: 'equilibrium', icon: 'arrowEquilibrium', title: 'Equilibrium arrow', tool: { name: 'arrow', style: 'equilibrium' } },
      { id: 'resonance', icon: 'arrowResonance', title: 'Resonance arrow', tool: { name: 'arrow', style: 'resonance' } },
      { id: 'retro', icon: 'arrowRetro', title: 'Retrosynthesis arrow', tool: { name: 'arrow', style: 'retro' } },
      { id: 'curved', icon: 'arrowCurved', title: 'Curved arrow (electron pair). Click an arrow to flip its curvature', tool: { name: 'arrow', style: 'curved' } },
      { id: 'fishhook', icon: 'arrowFishhook', title: 'Fishhook arrow (single electron)', tool: { name: 'arrow', style: 'fishhook' } },
      { id: 'noreaction', icon: 'arrowNoreaction', title: 'No reaction', tool: { name: 'arrow', style: 'noreaction' } },
      { id: 'dashedArrow', icon: 'arrowDashed', title: 'Dashed arrow', tool: { name: 'arrow', style: 'dashed' } },
    ],
  },
];

export function iconFor(def) {
  return def.iconHTML || ICONS[def.icon] || '';
}

/** Find the palette entry matching an editor tool object. */
export function matchToolId(tool) {
  for (const g of TOOL_GROUPS) {
    for (const t of g.tools) {
      const a = t.tool;
      if (a.name !== tool.name) continue;
      if (a.name === 'bond' && (a.order !== tool.order || a.stereo !== tool.stereo)) continue;
      if (a.name === 'ring' && a.template?.id !== tool.template?.id) continue;
      if (a.name === 'arrow' && a.style !== tool.style) continue;
      if (a.name === 'charge' && a.delta !== tool.delta) continue;
      return t.id;
    }
  }
  return null;
}
