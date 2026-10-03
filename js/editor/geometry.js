// Small 2D vector helpers (plain {x, y} objects).

export const add = (a, b) => ({ x: a.x + b.x, y: a.y + b.y });
export const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
export const mul = (a, k) => ({ x: a.x * k, y: a.y * k });
export const len = (a) => Math.hypot(a.x, a.y);
export const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
export const norm = (a) => { const l = len(a) || 1; return { x: a.x / l, y: a.y / l }; };
export const perp = (a) => ({ x: -a.y, y: a.x });
export const dot = (a, b) => a.x * b.x + a.y * b.y;
export const cross = (a, b) => a.x * b.y - a.y * b.x;
export const mid = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
export const angleOf = (v) => Math.atan2(v.y, v.x);
export const fromAngle = (t, l = 1) => ({ x: Math.cos(t) * l, y: Math.sin(t) * l });
export const deg = (r) => (r * 180) / Math.PI;
export const rad = (d) => (d * Math.PI) / 180;

export function rotateAround(p, c, t) {
  const cos = Math.cos(t), sin = Math.sin(t);
  const dx = p.x - c.x, dy = p.y - c.y;
  return { x: c.x + dx * cos - dy * sin, y: c.y + dx * sin + dy * cos };
}

/** Snap an angle (radians) to the nearest multiple of `stepDeg` degrees. */
export function snapAngle(t, stepDeg = 15) {
  const step = rad(stepDeg);
  return Math.round(t / step) * step;
}

/** Distance from point p to segment ab. */
export function distToSegment(p, a, b) {
  const ab = sub(b, a);
  const l2 = dot(ab, ab);
  if (l2 === 0) return dist(p, a);
  let t = dot(sub(p, a), ab) / l2;
  t = Math.max(0, Math.min(1, t));
  return dist(p, add(a, mul(ab, t)));
}

/** Normalise angle into (-PI, PI]. */
export function normAngle(t) {
  while (t <= -Math.PI) t += 2 * Math.PI;
  while (t > Math.PI) t -= 2 * Math.PI;
  return t;
}

export function pointInRect(p, r) {
  return p.x >= r.minX && p.x <= r.maxX && p.y >= r.minY && p.y <= r.maxY;
}

export function rectFromPoints(a, b) {
  return { minX: Math.min(a.x, b.x), minY: Math.min(a.y, b.y), maxX: Math.max(a.x, b.x), maxY: Math.max(a.y, b.y) };
}

/** Ray-casting point-in-polygon test. */
export function pointInPolygon(p, poly) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const a = poly[i], b = poly[j];
    if ((a.y > p.y) !== (b.y > p.y) && p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y) + a.x) inside = !inside;
  }
  return inside;
}
