import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Molecule } from '../js/core/molecule.js';
import {
  predictH1, predictC13, simulateSpectrum, symmetryClasses, formatH1Report, formatC13Report,
} from '../js/core/nmr.js';

// ---------------------------------------------------------------- helpers

/** Build a molecule from element symbols and [i, j, order] bonds (0-based indices). */
function build(els, bonds, coords = null) {
  const m = new Molecule();
  const ids = els.map((el, i) => m.addAtom({ el, x: coords?.[i]?.[0] ?? 0, y: coords?.[i]?.[1] ?? 0 }).id);
  for (const [a, b, o = 1] of bonds) m.addBond(ids[a], ids[b], o);
  return m;
}

/** Kekulé six-ring bonds over atom indices o..o+5. */
const kekule = (o) => [[o, o + 1, 2], [o + 1, o + 2, 1], [o + 2, o + 3, 2], [o + 3, o + 4, 1], [o + 4, o + 5, 2], [o + 5, o, 1]];
const C6 = ['C', 'C', 'C', 'C', 'C', 'C'];

const near = (actual, expected, tol, msg = '') => {
  assert.ok(Math.abs(actual - expected) <= tol, `${msg} expected ${expected}±${tol}, got ${actual}`);
};
/** Find the H signal whose shift lies in [lo, hi] (and optionally has nH). */
const sig = (res, lo, hi, nH = null) => {
  const s = res.signals.find((x) => x.shift >= lo && x.shift <= hi && (nH === null || x.nH === nH));
  assert.ok(s, `no signal in ${lo}-${hi}${nH ? ` with ${nH}H` : ''}: ${JSON.stringify(res.signals.map((x) => [x.shift, x.nH]))}`);
  return s;
};
const totalH = (res) => res.signals.reduce((s, x) => s + x.nH, 0);

const ethanol = () => build(['C', 'C', 'O'], [[0, 1], [1, 2]]);
const ethylAcetate = () => build(['C', 'C', 'O', 'O', 'C', 'C'], [[0, 1], [1, 2, 2], [1, 3], [3, 4], [4, 5]]);
const toluene = () => build([...C6, 'C'], [...kekule(0), [0, 6]]);
const benzene = () => build(C6, kekule(0));
const acetone = () => build(['C', 'C', 'O', 'C'], [[0, 1], [1, 2, 2], [1, 3]]);

// ---------------------------------------------------------------- 1H

test('ethanol: CH3 t, CH2 q, OH br s', () => {
  const r = predictH1(ethanol());
  assert.equal(r.signals.length, 3);
  assert.equal(totalH(r), 6);
  const ch3 = sig(r, 0.9, 1.5, 3);
  near(ch3.shift, 1.22, 0.3, 'CH3');
  assert.equal(ch3.multiplicity, 't');
  near(ch3.J[0], 7, 0.5);
  const ch2 = sig(r, 3.3, 4.0, 2);
  near(ch2.shift, 3.69, 0.3, 'CH2');
  assert.equal(ch2.multiplicity, 'q');
  const oh = r.signals.find((s) => s.exchangeable);
  assert.ok(oh);
  assert.equal(oh.nH, 1);
  assert.equal(oh.multiplicity, 'br s');
});

test('explicit hydrogens give the same result as implicit ones', () => {
  const m = build(['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'],
    [[0, 1], [1, 2], [0, 3], [0, 4], [0, 5], [1, 6], [1, 7], [2, 8]]);
  const a = predictH1(m).signals.map((s) => [s.shift, s.nH, s.multiplicity]);
  const b = predictH1(ethanol()).signals.map((s) => [s.shift, s.nH, s.multiplicity]);
  assert.deepEqual(a, b);
});

test('ethyl acetate: CH3CO s, OCH2 q, CH3 t', () => {
  const r = predictH1(ethylAcetate());
  assert.equal(r.signals.length, 3);
  const ac = sig(r, 1.8, 2.3, 3);
  assert.equal(ac.multiplicity, 's');
  near(ac.shift, 2.04, 0.3);
  const och2 = sig(r, 3.8, 4.4, 2);
  assert.equal(och2.multiplicity, 'q');
  near(och2.shift, 4.12, 0.3);
  const ch3 = sig(r, 1.0, 1.5, 3);
  assert.equal(ch3.multiplicity, 't');
  near(ch3.shift, 1.26, 0.3);
});

test('toluene: aromatic 7.0-7.35 (5H), CH3 s ~2.3', () => {
  const r = predictH1(toluene());
  const arom = r.signals.filter((s) => s.shift > 6.5);
  assert.equal(arom.reduce((s, x) => s + x.nH, 0), 5);
  assert.equal(arom.length, 3); // ortho, meta, para
  for (const s of arom) assert.ok(s.shift >= 7.0 && s.shift <= 7.35, `aromatic ${s.shift}`);
  const me = sig(r, 2.0, 2.6, 3);
  assert.equal(me.multiplicity, 's');
  near(me.shift, 2.34, 0.3);
});

test('benzene: single 6H singlet at 7.26', () => {
  const r = predictH1(benzene());
  assert.equal(r.signals.length, 1);
  assert.equal(r.signals[0].shift, 7.26);
  assert.equal(r.signals[0].nH, 6);
  assert.equal(r.signals[0].multiplicity, 's');
});

test('acetone: one 6H singlet ~2.17', () => {
  const r = predictH1(acetone());
  assert.equal(r.signals.length, 1);
  assert.equal(r.signals[0].nH, 6);
  assert.equal(r.signals[0].multiplicity, 's');
  assert.ok(r.signals[0].shift >= 2.0 && r.signals[0].shift <= 2.3);
});

test('2-chloropropane: CH septet ~4.1, CH3 doublet ~1.5', () => {
  const r = predictH1(build(['C', 'C', 'C', 'Cl'], [[0, 1], [1, 2], [1, 3]]));
  assert.equal(r.signals.length, 2);
  const ch = sig(r, 3.8, 4.5, 1);
  assert.equal(ch.multiplicity, 'sept');
  assert.equal(ch.peaks.length, 7);
  const me = sig(r, 1.2, 1.8, 6);
  assert.equal(me.multiplicity, 'd');
  near(me.shift, 1.55, 0.3);
});

test('acetaldehyde: CHO q ~9.8, CH3 d ~2.2', () => {
  const r = predictH1(build(['C', 'C', 'O'], [[0, 1], [1, 2, 2]]));
  const cho = sig(r, 9.4, 10.1, 1);
  assert.equal(cho.multiplicity, 'q');
  const me = sig(r, 1.9, 2.5, 3);
  assert.equal(me.multiplicity, 'd');
  near(me.J[0], cho.J[0], 1e-9);
});

test('styrene: three vinyl signals 5.2-6.7 with dd patterns', () => {
  const r = predictH1(build([...C6, 'C', 'C'], [...kekule(0), [0, 6], [6, 7, 2]]));
  const vinyl = r.signals.filter((s) => s.shift > 4.5 && s.shift < 6.9);
  assert.equal(vinyl.length, 3);
  for (const s of vinyl) {
    assert.ok(s.shift >= 5.1 && s.shift <= 6.8, `vinyl ${s.shift}`);
    assert.equal(s.nH, 1);
    assert.equal(s.multiplicity, 'dd');
  }
  const ch = vinyl[0];
  near(ch.shift, 6.72, 0.3);
  assert.deepEqual(ch.J, [17, 10]);
  assert.equal(totalH(r), 8);
});

test('alkene geometry from coordinates: E vs Z crotonaldehyde', () => {
  const els = ['C', 'C', 'C', 'C', 'O'];
  const bonds = [[0, 1], [1, 2, 2], [2, 3], [3, 4, 2]];
  const E = build(els, bonds, [[0, 0], [0.87, -0.5], [1.73, 0], [2.6, -0.5], [3.46, 0]]);
  const Z = build(els, bonds, [[0, -1], [0.87, -0.5], [1.73, 0], [1.73, 1], [2.6, 1.5]]);
  const jE = sig(predictH1(E), 5.8, 6.3, 1).J;
  const jZ = sig(predictH1(Z), 5.8, 6.3, 1).J;
  assert.ok(jE.includes(17), `E: ${jE}`);
  assert.ok(jZ.includes(10), `Z: ${jZ}`);
});

test('pyridine: α-H ~8.6, β ~7.25, γ ~7.65', () => {
  const r = predictH1(build(['N', 'C', 'C', 'C', 'C', 'C'], kekule(0)));
  assert.equal(r.signals.length, 3);
  near(sig(r, 8.3, 8.9, 2).shift, 8.6, 0.15);
  near(sig(r, 7.5, 7.8, 1).shift, 7.65, 0.15);
  near(sig(r, 7.0, 7.45, 2).shift, 7.25, 0.15);
});

test('benzaldehyde: CHO ~10.0', () => {
  const r = predictH1(build([...C6, 'C', 'O'], [...kekule(0), [0, 6], [6, 7, 2]]));
  near(r.signals[0].shift, 10.0, 0.2);
  assert.equal(r.signals[0].nH, 1);
  near(sig(r, 7.6, 8.1, 2).shift, 7.87, 0.3); // ortho H
});

test('acetic acid: COOH > 10 (br s, exchangeable), CH3 ~2.1', () => {
  const r = predictH1(build(['C', 'C', 'O', 'O'], [[0, 1], [1, 2, 2], [1, 3]]));
  const oh = r.signals[0];
  assert.ok(oh.shift > 10);
  assert.equal(oh.exchangeable, true);
  assert.equal(oh.multiplicity, 'br s');
  near(sig(r, 1.8, 2.4, 3).shift, 2.1, 0.2);
});

test('exchangeable NH/OH do not split neighbours (methylamine, phenol)', () => {
  const r = predictH1(build(['C', 'N'], [[0, 1]]));
  const me = sig(r, 2.0, 2.8, 3);
  assert.equal(me.multiplicity, 's');
  const nh = r.signals.find((s) => s.exchangeable);
  assert.equal(nh.nH, 2);
  const ph = predictH1(build([...C6, 'O'], [...kekule(0), [0, 6]]));
  const phOH = ph.signals.find((s) => s.exchangeable);
  assert.ok(phOH.shift >= 4.5 && phOH.shift <= 7);
});

test('1-propanol: middle CH2 is a sextet', () => {
  const r = predictH1(build(['C', 'C', 'C', 'O'], [[0, 1], [1, 2], [2, 3]]));
  assert.equal(sig(r, 1.4, 1.9, 2).multiplicity, 'sext');
  assert.equal(sig(r, 3.4, 3.9, 2).multiplicity, 't');
});

test('peaks: Pascal intensities sum to nH and spacing = J / frequency', () => {
  const r = predictH1(ethanol(), { frequency: 500 });
  const ch2 = sig(r, 3.3, 4.0, 2);
  assert.equal(ch2.peaks.length, 4);
  const sum = ch2.peaks.reduce((s, p) => s + p.intensity, 0);
  near(sum, 2, 1e-9);
  const ratios = ch2.peaks.map((p) => p.intensity / ch2.peaks[0].intensity);
  ratios.forEach((v, i) => near(v, [1, 3, 3, 1][i], 1e-9));
  near(ch2.peaks[1].shift - ch2.peaks[0].shift, ch2.J[0] / 500, 1e-4);
});

// ---------------------------------------------------------------- 13C

test('13C benzene: one signal, count 6, 128.5', () => {
  const r = predictC13(benzene());
  assert.equal(r.signals.length, 1);
  assert.equal(r.signals[0].count, 6);
  assert.equal(r.signals[0].shift, 128.5);
  assert.equal(r.signals[0].multiplicity, 'CH');
});

test('13C toluene: 5 signals incl. ~21 (CH3) and ~137 (ipso)', () => {
  const r = predictC13(toluene());
  assert.equal(r.signals.length, 5);
  const me = r.signals.find((s) => s.multiplicity === 'CH3');
  near(me.shift, 21.4, 5);
  const ipso = r.signals.find((s) => s.multiplicity === 'C');
  near(ipso.shift, 137.8, 5);
  assert.deepEqual(r.signals.map((s) => s.count).sort(), [1, 1, 1, 2, 2]);
});

test('13C acetone: C=O ~206, CH3 ~30', () => {
  const r = predictC13(acetone());
  assert.equal(r.signals.length, 2);
  near(r.signals[0].shift, 206.7, 5);
  near(r.signals[1].shift, 30.8, 5);
  assert.equal(r.signals[1].count, 2);
});

test('13C ethanol and ethyl acetate within 5 ppm of literature', () => {
  const e = predictC13(ethanol()).signals.map((s) => s.shift);
  near(e[0], 58.3, 5);
  near(e[1], 18.4, 5);
  const ea = predictC13(ethylAcetate()).signals.map((s) => s.shift);
  [171.1, 60.4, 21.0, 14.2].forEach((v, i) => near(ea[i], v, 5, `EtOAc C${i}`));
});

test('13C pyridine, styrene and acetonitrile', () => {
  const py = predictC13(build(['N', 'C', 'C', 'C', 'C', 'C'], kekule(0))).signals.map((s) => s.shift);
  [149.9, 135.9, 123.8].forEach((v, i) => near(py[i], v, 3));
  const st = predictC13(build([...C6, 'C', 'C'], [...kekule(0), [0, 6], [6, 7, 2]]));
  const ch2 = st.signals.find((s) => s.multiplicity === 'CH2');
  near(ch2.shift, 113.7, 5);
  const mecn = predictC13(build(['C', 'C', 'N'], [[0, 1], [1, 2, 3]])).signals;
  near(mecn[0].shift, 117.7, 3);
});

// ---------------------------------------------------------------- utilities

test('symmetryClasses: isopropyl methyls and benzene ortho/meta pairs are equivalent', () => {
  const ipa = build(['C', 'C', 'C', 'O'], [[0, 1], [1, 2], [1, 3]]);
  const ids = [...ipa.atoms.keys()];
  const cls = symmetryClasses(ipa);
  assert.equal(cls.get(ids[0]), cls.get(ids[2]));
  assert.notEqual(cls.get(ids[0]), cls.get(ids[1]));
  const t = toluene();
  const tid = [...t.atoms.keys()];
  const tc = symmetryClasses(t);
  assert.equal(tc.get(tid[1]), tc.get(tid[5])); // ortho
  assert.equal(tc.get(tid[2]), tc.get(tid[4])); // meta
  assert.notEqual(tc.get(tid[1]), tc.get(tid[2]));
  assert.equal(new Set(tc.values()).size, 5);
});

test('simulateSpectrum: normalised Lorentzian on requested axis', () => {
  const { x, y } = simulateSpectrum([{ shift: 1, intensity: 3 }, { shift: 2, intensity: 1 }], { min: 0, max: 3, points: 3001 });
  assert.ok(x instanceof Float64Array && y instanceof Float64Array);
  assert.equal(x.length, 3001);
  near(x[0], 0, 1e-12);
  near(x[3000], 3, 1e-12);
  near(Math.max(...y), 1, 1e-12);
  near(y[1000], 1, 1e-9); // tallest line at 1.00 ppm
  near(y[2000], 1 / 3, 0.01);
});

test('report formatting', () => {
  const h = formatH1Report(predictH1(ethylAcetate()));
  assert.equal(h, '1H NMR (400 MHz, CDCl3) δ 4.12 (q, J = 7.0 Hz, 2H), 2.03 (s, 3H), 1.26 (t, J = 7.0 Hz, 3H).');
  const c = formatC13Report(predictC13(acetone()), { solvent: 'DMSO-d6' });
  assert.match(c, /^13C NMR \(100 MHz, DMSO-d6\) δ 20\d\.\d, \d+\.\d\.$/);
  const fake = { signals: [
    { shift: 7.30, nH: 2, multiplicity: 'm', J: [] },
    { shift: 7.20, nH: 3, multiplicity: 'm', J: [] },
    { shift: 2.35, nH: 3, multiplicity: 's', J: [] },
  ] };
  assert.equal(formatH1Report(fake, { frequency: 500 }), '1H NMR (500 MHz, CDCl3) δ 7.30–7.20 (m, 5H), 2.35 (s, 3H).');
});
