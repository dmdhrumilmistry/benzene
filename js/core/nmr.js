// Empirical 1H / 13C NMR shift prediction.
//
// This is an additive-increment estimator in the spirit of ChemDraw's
// "ChemNMR basic" mode: every H / C is classified by its local environment and
// a base value is corrected by tabulated substituent increments
// (Shoolery / Pascual–Meier–Simon for 1H, Grant–Paul and benzene / alkene
// increment tables for 13C). Numbers are taken (and rounded) from the
// standard compilations (Pretsch, "Structure Determination of Organic
// Compounds"; Silverstein). Expect ±0.3 ppm (1H) / ±5 ppm (13C) for ordinary
// organic molecules and larger errors for exotic ones.
//
// Contributors: all increments live in the tables at the top of this file.
// Every table is keyed by a *substituent group key* produced by
// `Ctx.group(fromAtom, substituentAtom)`; see GROUP KEYS below. If a key is
// missing from a table the FALLBACK chain is followed (e.g. 'NHR' -> 'NR2').

import { analyzeRings } from './rings.js';

// ===========================================================================
// GROUP KEYS (substituent classes, as seen from the atom they are attached to)
//
//   H, F, Cl, Br, I
//   CH3 | alkyl (CH2R) | alkylCH (CHR2) | alkylC (CR3)    plain sp3 carbon
//   alkylO | alkylN | alkylHal | CX3                     sp3 carbon bearing O / N / halogen / >=2 halogens
//   aryl | vinyl | alkynyl
//   fused (aromatic ring-fusion neighbour; fused5 / fused5X when the fused ring is
//   a pyrrole / furan ring and the neighbour is its C / heteroatom, as in indole)
//   CN (nitrile) | CHO | COR (ketone) | COOH | COO- | COOR | CONR2 | COCl | C=N
//   NH2 | NHR | NR2 (amines) | NAr (N-aryl amine) | Namide (N-C=O / N-SO2)
//   NR3+ | NO2 | Narom (aromatic N, e.g. pyrrole N-R) | N=X (imine / azo N)
//   OH | O- | OR | Ovinyl | OAr | OCOR (ester O) | OCOAr (ester O of an aryl ester) | OSO2
//   SH | SR | SOR | SO2R | SiR3 | P | B | =O (exocyclic double-bonded O/S) | other
// ===========================================================================

const FALLBACK = {
  CH3: 'alkyl', alkylCH: 'alkyl', alkylC: 'alkyl',
  alkylO: 'alkyl', alkylN: 'alkyl', alkylHal: 'alkyl', CX3: 'alkylHal',
  fused: 'aryl', fused5: 'fused', fused5X: 'fused', OCOAr: 'OCOR', 'COO-': 'COOH', COCl: 'COR', 'C=N': 'COR',
  NH2: 'NR2', NHR: 'NR2', NAr: 'NR2', 'NR3+': 'NR2', Narom: 'NAr', 'N=X': 'NAr',
  'O-': 'OH', Ovinyl: 'OR', OSO2: 'OCOR', SH: 'SR', SOR: 'SO2R',
};

/** Plain alkyl-type groups (contribute nothing as an alpha group to sp3 1H). */
const ALKYL = new Set(['CH3', 'alkyl', 'alkylCH', 'alkylC', 'alkylO', 'alkylN', 'alkylHal', 'CX3']);
const PLAIN_ALKYL = new Set(['CH3', 'alkyl', 'alkylCH', 'alkylC']);
const HALOGENS = new Set(['F', 'Cl', 'Br', 'I']);
const CONJ = new Set(['aryl', 'fused', 'vinyl']);

// ---------------------------------------------------------------------------
// 1H: sp3 C–H.  δ = base(CHn) + Σ α-increments + Σ β-increments
// Base values are for CH3/CH2/CH flanked only by alkyl groups.
const H_SP3_BASE = { 4: 0.23, 3: 0.86, 2: 1.37, 1: 1.50 };

// α-increments [CH3, CH2, CH]: δ(H3C–X) − 0.86, δ(CH3CH2–X) − 1.37, δ((CH3)2CH–X) − 1.50.
const H_SP3_ALPHA = {
  alkyl: [0, 0, 0],
  vinyl: [0.85, 0.67, 0.80],
  aryl: [1.49, 1.26, 1.40],
  alkynyl: [0.94, 0.79, 0.90],
  CHO: [1.34, 1.09, 0.92],
  COR: [1.31, 1.07, 1.10],
  COOH: [1.24, 1.01, 1.07],
  'COO-': [1.05, 0.85, 0.90],
  COOR: [1.17, 0.95, 1.00],
  CONR2: [1.16, 0.86, 0.90],
  COCl: [1.80, 1.50, 1.40],
  'C=N': [1.00, 0.90, 0.90],
  CN: [1.14, 0.99, 1.20],
  NH2: [1.61, 1.37, 1.60],
  NHR: [1.57, 1.28, 1.40],
  NR2: [1.36, 1.16, 1.30],
  NAr: [2.00, 1.75, 2.10],
  Namide: [1.94, 1.88, 2.60],
  'NR3+': [2.44, 2.00, 2.00],
  NO2: [3.43, 3.00, 3.17],
  Narom: [2.75, 2.60, 3.00],
  'N=X': [2.40, 2.20, 2.30],
  OH: [2.53, 2.32, 2.51],
  OR: [2.40, 2.11, 2.15],
  OAr: [2.94, 2.65, 3.00],
  OCOR: [2.81, 2.75, 3.50],
  OCOAr: [3.05, 2.95, 3.75],
  OSO2: [2.85, 2.70, 3.30],
  F: [3.40, 3.18, 3.34],
  Cl: [2.19, 2.20, 2.64],
  Br: [1.82, 2.06, 2.71],
  I: [1.30, 1.83, 2.74],
  SR: [1.23, 1.19, 1.50],
  SOR: [1.76, 1.30, 1.30],
  SO2R: [2.14, 1.60, 1.60],
  SiR3: [-0.86, -0.85, -0.80],
  P: [0.90, 0.80, 0.80],
  B: [-0.60, -0.50, -0.50],
};

// β-increments (group on the neighbouring sp3 carbon), all CHn types.
const H_SP3_BETA = {
  vinyl: 0.10, aryl: 0.35, alkynyl: 0.25, CHO: 0.25, COR: 0.20, COOH: 0.30, COOR: 0.25,
  CONR2: 0.20, COCl: 0.40, CN: 0.45, NH2: 0.20, NR2: 0.15, NAr: 0.30, Namide: 0.25,
  'NR3+': 0.50, NO2: 0.70, Narom: 0.50, 'N=X': 0.30, OH: 0.30, OR: 0.25, OAr: 0.45,
  OCOR: 0.40, OSO2: 0.45, F: 0.45, Cl: 0.55, Br: 0.80, I: 1.00, SR: 0.40, SOR: 0.40,
  SO2R: 0.45, SiR3: -0.30,
};

// ---------------------------------------------------------------------------
// 1H: benzene ring.  δ = 7.26 + Σ [ortho, meta, para] increments.
const AROMATIC_H = {
  CH3: [-0.18, -0.11, -0.21],
  alkyl: [-0.14, -0.06, -0.17],
  alkylCH: [-0.13, -0.08, -0.18],
  alkylC: [0.02, -0.08, -0.21],
  alkylO: [-0.07, -0.07, -0.07],
  alkylN: [-0.05, -0.05, -0.05],
  alkylHal: [0.00, 0.00, 0.00],
  CX3: [0.32, 0.14, 0.20],
  fused5: [0.30, -0.10, -0.10], // benzo ring of indole / benzofuran, fusion C3a side
  fused5X: [0.10, 0.10, -0.15], // ... fusion C7a side (ring N / O)
  vinyl: [0.06, -0.03, -0.10],
  alkynyl: [0.15, -0.02, -0.01],
  aryl: [0.37, 0.20, 0.10],
  F: [-0.26, 0.00, -0.20],
  Cl: [0.03, -0.02, -0.09],
  Br: [0.18, -0.08, -0.04],
  I: [0.39, -0.21, 0.00],
  OH: [-0.56, -0.12, -0.45],
  OR: [-0.48, -0.09, -0.44],
  OAr: [-0.29, -0.05, -0.23],
  OCOR: [-0.25, 0.03, -0.13],
  OSO2: [-0.05, 0.07, -0.01],
  NH2: [-0.75, -0.25, -0.65],
  NHR: [-0.72, -0.20, -0.66],
  NR2: [-0.66, -0.18, -0.67],
  NAr: [-0.25, 0.00, -0.35],
  Namide: [0.12, -0.07, -0.28],
  'NR3+': [0.40, 0.20, 0.30],
  NO2: [0.95, 0.26, 0.38],
  Narom: [0.15, 0.15, 0.00],
  'N=X': [0.67, 0.20, 0.20],
  CN: [0.36, 0.18, 0.28],
  CHO: [0.56, 0.22, 0.29],
  COR: [0.62, 0.14, 0.21],
  COOH: [0.85, 0.18, 0.27],
  COOR: [0.71, 0.10, 0.21],
  CONR2: [0.61, 0.10, 0.17],
  COCl: [0.84, 0.22, 0.36],
  'C=N': [0.50, 0.15, 0.15],
  SH: [-0.08, -0.16, -0.22],
  SR: [-0.08, -0.10, -0.24],
  SOR: [0.38, 0.20, 0.26],
  SO2R: [0.60, 0.25, 0.33],
  SiR3: [0.22, -0.02, -0.02],
  P: [0.30, 0.10, 0.10],
  B: [0.50, 0.10, 0.15],
};

// Ring nitrogen effects in six-membered heteroaromatics (pyridine: α 8.60, β 7.25, γ 7.64).
// [α, β, γ]; a second α-nitrogen counts at ALPHA_N_DAMP of the first.
const PYRIDINE_N_H = [1.34, 0.00, 0.38];
const PYRIDINE_N_C = [21.4, -4.7, 7.4];
const ALPHA_N_DAMP_H = 0.5;
const ALPHA_N_DAMP_C = 0.4;

// Five-membered heteroaromatics: base shift by the π-donor heteroatom and
// position (a = adjacent to it, b = one further).
const FIVE_RING_H = {
  N: { a: 6.68, b: 6.22 }, // pyrrole
  O: { a: 7.42, b: 6.38 }, // furan
  S: { a: 7.31, b: 7.10 }, // thiophene
  Se: { a: 7.90, b: 7.25 },
};
const FIVE_RING_C = {
  N: { a: 118.2, b: 108.0 },
  O: { a: 142.7, b: 109.6 },
  S: { a: 125.4, b: 127.2 },
  Se: { a: 129.8, b: 131.0 },
};
// Extra pyridine-type (=N–) nitrogen in a five-ring (imidazole, oxazole, thiazole…): [adjacent, 1,3].
const AZOLE_N_H = [1.10, 0.30];
const AZOLE_N_C = [17.0, 4.0];

// ---------------------------------------------------------------------------
// 1H: alkenes (Pascual–Meier–Simon).  δ = 5.25 + Z_gem + Z_cis + Z_trans.
const VINYL_H = {
  alkyl: [0.45, -0.22, -0.28],
  alkylO: [0.64, -0.01, -0.02],
  alkylN: [0.58, -0.10, -0.08],
  alkylHal: [0.70, 0.11, -0.04],
  CX3: [0.66, 0.61, 0.32],
  vinyl: [1.00, -0.09, -0.23],
  alkynyl: [0.47, 0.38, 0.12],
  aryl: [1.38, 0.36, -0.07],
  CHO: [1.02, 0.95, 1.17],
  COR: [1.10, 1.12, 0.87],
  COOH: [0.97, 1.41, 0.71],
  COOR: [0.80, 1.18, 0.55],
  CONR2: [1.37, 0.98, 0.46],
  COCl: [1.11, 1.46, 1.01],
  CN: [0.27, 0.75, 0.55],
  F: [1.54, -0.40, -1.02],
  Cl: [1.08, 0.18, 0.13],
  Br: [1.07, 0.45, 0.55],
  I: [1.14, 0.81, 0.88],
  OR: [1.22, -1.07, -1.21],
  OAr: [1.21, -0.60, -1.00],
  OCOR: [2.11, -0.35, -0.64],
  NR2: [0.80, -1.26, -1.21],
  Namide: [2.08, -0.57, -0.72],
  Narom: [1.80, -0.60, -0.80],
  NO2: [1.87, 1.32, 0.62],
  SR: [1.11, -0.29, -0.13],
  SOR: [1.27, 0.67, 0.41],
  SO2R: [1.55, 1.16, 0.93],
  SiR3: [0.90, 0.90, 0.60],
  P: [0.66, 0.88, 0.67],
};

// 1H: terminal alkyne ≡C–H. δ = 1.90 + increment of the group on the other sp carbon.
const ALKYNE_H_BASE = 1.90;
const ALKYNE_H = {
  alkyl: 0, aryl: 1.15, vinyl: 1.00, alkynyl: 0.10, alkylO: 0.55, alkylN: 0.35, alkylHal: 0.50,
  CHO: 1.40, COR: 1.35, COOH: 1.10, COOR: 1.05, CONR2: 0.90, CN: 0.90, SiR3: 0.50,
};

// 1H coupling constants (Hz).
const J = {
  alkyl: 7.0, // H–C–C–H, freely rotating sp3
  allylic: 6.5, // H–C(sp2)–C(sp3)–H
  aldehydeSp3: 2.5, // HC(=O)–C(sp3)–H
  aldehydeSp2: 7.8, // HC(=O)–C(sp2)–H
  diene: 10.5, // H–C(sp2)–C(sp2)–H across a single bond
  cis: 10.0,
  trans: 17.0,
  gem: 1.5, // =CH2
  ortho: 8.0,
  orthoNextToN: 5.0, // pyridine J(2,3)
  fiveRing: 3.5,
};
const J_MERGE_TOL = 1.0; // couplings closer than this are treated as equal (n+1 rule)
const STRONG_COUPLING_RATIO = 4; // Δν/J below this => second order => 'm'
const MULT_NAMES = ['s', 'd', 't', 'q', 'quint', 'sext', 'sept'];

// ---------------------------------------------------------------------------
// 13C: sp3 carbon (Grant–Paul).
// δ = -2.3 + 9.1 nα + 9.4 nβ − 2.5 nγ + 0.3 nδ + branching + Σ group increments.
const C13_SP3 = { base: -2.3, alpha: 9.1, beta: 9.4, gamma: -2.5, delta: 0.3 };
// Branching corrections [observed carbon degree][neighbour degree]; degree = # carbon neighbours.
const C13_BRANCH = {
  1: { 3: -1.1, 4: -3.4 },
  2: { 3: -2.5, 4: -7.2 },
  3: { 2: -3.7, 3: -9.5, 4: -15.0 },
  4: { 1: -1.5, 2: -8.4, 3: -15.0, 4: -25.0 },
};
// Ring corrections for sp3 carbons in small / medium rings (by smallest ring size).
const C13_RING = { 3: -18.7, 4: -3.0, 5: -9.0, 6: -5.0, 7: -1.0 };
// Functional-group increments [α, β, γ] (group replaces an H on the α/β/γ carbon).
const C13_SP3_INC = {
  aryl: [23, 9, -2],
  vinyl: [20, 7, -2],
  alkynyl: [4.5, 5.5, -3.5],
  CN: [4, 3, -3],
  CHO: [31, 0, -2],
  COR: [31, 1, -2],
  COOH: [21, 3, -2],
  'COO-': [24, 3, -2],
  COOR: [22.5, 3, -2],
  CONR2: [25, 2.5, -0.5],
  COCl: [33, 2, -3.5],
  'C=N': [24, 2, -2],
  NH2: [29, 11, -5],
  NHR: [37, 8, -4],
  NR2: [42, 6, -3],
  NAr: [33, 8, -4],
  Namide: [30, 7, -5],
  'NR3+': [31, 5, -7],
  NO2: [63, 4, -4.6],
  Narom: [38, 8, -5],
  'N=X': [50, 7, -4],
  OH: [49, 10, -6],
  OR: [58, 8, -4],
  OAr: [57, 7, -4],
  OCOR: [54, 6, -6],
  OSO2: [58, 7, -6],
  F: [70, 8, -7],
  Cl: [31, 11, -4],
  Br: [20, 11, -3],
  I: [-6, 11, -1],
  SH: [11, 12, -4],
  SR: [20, 7, -3],
  SOR: [43, 0, -3],
  SO2R: [45, 0, -3],
  SiR3: [2, 3, -1],
  P: [10, 2, -2],
  B: [8, 2, 0],
};
// Extra α-effect of heavy halogens on branched carbons, per additional carbon neighbour.
const C13_HALOGEN_BRANCH = { Cl: 6, Br: 9, I: 11 };
const HETERO_DAMP = 0.8; // each further α (or β) heteroatom group counts 0.8× the previous

// 13C: benzene ring. δ = 128.5 + Σ [ipso, ortho, meta, para].
const AROMATIC_C = {
  CH3: [9.3, 0.7, -0.1, -2.9],
  alkyl: [15.6, -0.5, 0.0, -2.6],
  alkylCH: [20.1, -2.0, 0.0, -2.5],
  alkylC: [22.2, -3.4, -0.4, -3.1],
  alkylO: [12.4, -1.2, 0.2, -1.1],
  alkylN: [15.0, -1.4, -0.2, -2.0],
  alkylHal: [9.1, 0.0, 0.2, -0.2],
  CX3: [2.6, -3.1, 0.4, 3.4],
  vinyl: [8.9, -2.3, -0.1, -0.8],
  alkynyl: [-6.2, 3.6, -0.4, -0.3],
  aryl: [13.1, -1.1, 0.5, -1.1],
  fused: [5.5, -0.5, 0.0, -1.5],
  fused5: [16.8, -7.7, 0.0, -6.5], // fitted to indole
  fused5X: [15.0, -17.4, 0.0, -8.7],
  F: [34.8, -13.0, 1.6, -4.4],
  Cl: [6.3, 0.4, 1.4, -1.9],
  Br: [-5.8, 3.2, 1.6, -1.6],
  I: [-34.1, 8.9, 1.6, -1.1],
  OH: [26.9, -12.8, 1.4, -7.4],
  'O-': [39.6, -8.2, 1.9, -13.6],
  OR: [31.4, -14.4, 1.0, -7.7],
  OAr: [29.1, -9.5, 0.3, -5.3],
  OCOR: [22.4, -7.1, 0.4, -3.2],
  OSO2: [21.2, -6.1, 1.1, -1.5],
  NH2: [18.2, -13.4, 0.8, -10.0],
  NHR: [20.7, -16.2, 0.7, -11.8],
  NR2: [22.4, -15.7, 0.8, -11.8],
  NAr: [14.7, -10.6, 0.9, -7.6],
  Namide: [9.7, -8.1, 0.2, -4.4],
  'NR3+': [2.5, -6.0, 1.5, 1.5],
  NO2: [19.9, -4.9, 0.9, 6.1],
  Narom: [12.3, -7.9, 1.1, -2.8],
  'N=X': [24.1, -5.7, 0.5, 2.4],
  CN: [-15.7, 3.6, 0.7, 4.3],
  CHO: [8.2, 1.2, 0.5, 5.8],
  COR: [8.9, 0.0, -0.1, 4.4],
  COOH: [2.1, 1.6, -0.1, 5.2],
  'COO-': [8.0, 1.0, 0.0, 3.0],
  COOR: [2.0, 1.2, -0.1, 4.3],
  CONR2: [5.0, -1.2, 0.1, 3.4],
  COCl: [4.7, 2.7, 0.3, 6.6],
  'C=N': [7.5, -0.5, 0.0, 2.5],
  SH: [2.1, 0.7, 0.3, -3.2],
  SR: [10.0, -1.9, 0.2, -3.6],
  SOR: [17.1, -5.0, 0.8, 2.5],
  SO2R: [12.3, -1.4, 0.8, 5.1],
  SiR3: [11.6, 4.9, -0.7, 0.4],
  P: [8.8, 5.2, 0.0, 0.2],
  B: [1.0, 7.0, -0.5, 3.0],
};

// 13C: alkenes. δ = 123.3 + Σ [same carbon, other carbon] increments.
// Plain alkyl chains are counted instead: same carbon +10.6 (α) +7.2 (β) −1.5 (γ);
// other carbon −7.9 (α') −1.8 (β') +1.5 (γ').
const ALKENE_C_BASE = 123.3;
const ALKENE_C_ALKYL = { a: 10.6, b: 7.2, g: -1.5, a2: -7.9, b2: -1.8, g2: 1.5 };
const ALKENE_C = {
  alkylO: [14.2, -8.4],
  alkylN: [15.0, -10.5],
  alkylHal: [10.2, -6.0],
  CX3: [5.0, -1.0],
  aryl: [12.5, -11.0],
  vinyl: [13.6, -7.0],
  alkynyl: [-4.0, 6.0],
  CHO: [13.1, 12.7],
  COR: [15.0, 5.8],
  COOH: [4.2, 8.9],
  'COO-': [8.0, 4.0],
  COOR: [6.3, 7.0],
  CONR2: [8.0, 4.0],
  COCl: [8.1, 14.0],
  CN: [-15.1, 14.2],
  'C=N': [12.0, 0.0],
  F: [24.9, -34.3],
  Cl: [2.8, -6.1],
  Br: [-8.6, -0.9],
  I: [-38.1, 7.0],
  OR: [29.0, -39.0],
  OAr: [28.0, -33.0],
  OCOR: [18.4, -26.7],
  NR2: [28.0, -32.0],
  Namide: [6.2, -28.8],
  Narom: [6.0, -24.0],
  NO2: [22.3, -0.9],
  SR: [12.6, -13.7],
  SOR: [18.0, 0.0],
  SO2R: [14.3, 7.9],
  SiR3: [16.9, 6.7],
  P: [9.0, 9.0],
};

// 13C: alkynes. δ = 71.9 + Σ [same carbon, other carbon].
const ALKYNE_C_BASE = 71.9;
const ALKYNE_C_ALKYL = { a: 6.9, b: 4.8, a2: -5.7, b2: 2.3 };
const ALKYNE_C = {
  aryl: [11.8, 5.3],
  vinyl: [11.0, 6.0],
  alkynyl: [-4.0, 4.0],
  alkylO: [10.0, 2.0],
  alkylN: [11.0, 0.0],
  alkylHal: [8.0, 3.0],
  CX3: [0.0, 0.0],
  CHO: [10.0, 12.0],
  COR: [9.0, 10.0],
  COOH: [2.0, 7.0],
  COOR: [2.0, 4.0],
  CONR2: [5.0, 4.0],
  CN: [-15.0, 10.0],
  OR: [18.0, -50.0],
  SiR3: [18.0, 20.0],
};

// 13C: carbonyl and related sp2 carbons (base values; see c13Carbonyl()).
const C13_CO = {
  ketone: 206.5, ketoneBeta: 2.5, ketoneConj: -8.5, ketoneConj2: -2.0, ketoneAlkynyl: -22, ketoneDione: -8,
  aldehyde: 200.5, aldehydeBeta: 2.0, aldehydeConj: 192.0, aldehydeVinyl: 193.5, formaldehyde: 194.0,
  formate: 161.0, formamide: 162.5,
  acid: 177.5, acidConj: 172.0, carboxylate: 178.0,
  ester: 171.0, esterConj: 166.5, anhydride: 167.0, anhydrideConj: 163.0,
  amide: 170.5, amideConj: 167.5,
  acidHalide: 170.0, acidHalideConj: 168.0,
  thioester: 197.0, thioesterConj: 191.0,
  betaPerC: 2.0, // acid / ester / amide: per carbon on the α-carbon (max 3)
  ring5: 7.0, ring4: -3.0, // ring-strain corrections (cyclopentanone, β-lactam…)
  carbonate: 155.5, carbamate: 156.5, urea: 158.5,
  thioamide: 200.0, thiourea: 182.0, thione: 210.0,
  aldimine: 165.0, aldimineConj: 160.5, ketimine: 170.0, ketimineConj: 166.0,
  aldoxime: 149.0, ketoxime: 157.0, hydrazone: 143.0, amidine: 165.0, guanidine: 157.0, imidate: 163.0,
  nitrile: 118.0, nitrileConj: 118.8,
  alleneCentral: 210.0, alleneTerminal: 75.0, isocyanate: 123.0, co2: 125.0, ketene: 194.0,
};

// Readable labels for environment strings.
const GROUP_LABEL = {
  CH3: 'CH3', alkyl: 'alkyl', alkylCH: 'alkyl', alkylC: 'alkyl', alkylO: 'CH2O', alkylN: 'CH2N',
  alkylHal: 'CH2X', CX3: 'CX3', aryl: 'aryl', fused: 'fused ring', fused5: 'fused ring', fused5X: 'fused ring', OCOAr: 'OC(O)Ar', vinyl: 'C=C', alkynyl: 'C≡C',
  CN: 'C≡N', CHO: 'CHO', COR: 'C=O', COOH: 'COOH', 'COO-': 'CO2−', COOR: 'CO2R', CONR2: 'C(O)N',
  COCl: 'C(O)X', 'C=N': 'C=N', NH2: 'NH2', NHR: 'NHR', NR2: 'NR2', NAr: 'N-aryl', Namide: 'N-C(O)',
  'NR3+': 'N+', NO2: 'NO2', Narom: 'N (aromatic)', 'N=X': 'N=', OH: 'OH', 'O-': 'O−', OR: 'OR',
  Ovinyl: 'O-vinyl', OAr: 'OAr', OCOR: 'OC(O)R', OSO2: 'OSO2R', SH: 'SH', SR: 'SR', SOR: 'S(O)R',
  SO2R: 'SO2R', SiR3: 'SiR3', P: 'P', B: 'B', F: 'F', Cl: 'Cl', Br: 'Br', I: 'I',
};
const POS_NAME = ['ipso', 'ortho', 'meta', 'para'];

// ===========================================================================
// Helpers
// ===========================================================================

function lookup(table, key) {
  let k = key;
  for (let i = 0; i < 8 && k; i++) {
    if (Object.prototype.hasOwnProperty.call(table, k)) return table[k];
    k = FALLBACK[k];
  }
  return null;
}

/** Sum with diminishing weight for further (positive) contributions. */
function damped(values, f = HETERO_DAMP) {
  const pos = values.filter((v) => v > 0).sort((a, b) => b - a);
  let s = 0;
  pos.forEach((v, i) => { s += v * f ** i; });
  for (const v of values) if (v <= 0) s += v;
  return s;
}

const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
const round = (v, d) => Math.round(v * 10 ** d) / 10 ** d;
const label = (g) => GROUP_LABEL[g] ?? g;

/** Topological symmetry classes via iterative (Morgan-style) refinement. */
function computeSymmetry(mol, aromBonds, atomRings) {
  const isSuppressedH = (a) => a.el === 'H' && mol.neighbors(a.id).some((n) => mol.getAtom(n).el !== 'H');
  const nodes = mol.atomList().filter((a) => !isSuppressedH(a)).map((a) => a.id);
  const nodeSet = new Set(nodes);
  const bcode = (b) => (aromBonds.has(b.id) ? 'a' : String(b.order));
  const nbrs = new Map(nodes.map((id) => [id, mol.bondsOf(id)
    .map((b) => [mol.otherAtom(b, id), bcode(b)]).filter(([n]) => nodeSet.has(n))]));
  const rank = (inv) => {
    const keys = [...new Set(inv.values())].sort();
    const idx = new Map(keys.map((k, i) => [k, i]));
    const out = new Map();
    for (const [id, k] of inv) out.set(id, idx.get(k));
    return { cls: out, count: keys.length };
  };
  const init = new Map();
  for (const id of nodes) {
    const a = mol.getAtom(id);
    const bonds = nbrs.get(id).map(([, c]) => c).sort().join('');
    init.set(id, [a.el, a.isotope ?? '', a.charge, a.radical, nbrs.get(id).length, mol.totalH(id),
      aromBonds.size && mol.bondsOf(id).some((b) => aromBonds.has(b.id)) ? 1 : 0,
      atomRings.has(id) ? 1 : 0, bonds].join('|'));
  }
  let { cls, count } = rank(init);
  for (let iter = 0; iter < nodes.length; iter++) {
    const inv = new Map();
    for (const id of nodes) {
      const env = nbrs.get(id).map(([n, c]) => `${c}${cls.get(n)}`).sort().join(',');
      inv.set(id, `${cls.get(id)}:${env}`);
    }
    const next = rank(inv);
    cls = next.cls;
    if (next.count === count) break;
    count = next.count;
  }
  // Suppressed (explicit, attached) hydrogens get a class derived from their parent.
  const result = new Map(cls);
  for (const a of mol.atoms.values()) {
    if (!isSuppressedH(a)) continue;
    const parent = mol.neighbors(a.id).find((n) => mol.getAtom(n).el !== 'H');
    result.set(a.id, count + cls.get(parent));
  }
  return result;
}

/** Per-molecule perception cache + group classification. */
class Ctx {
  constructor(mol) {
    this.mol = mol;
    const r = analyzeRings(mol);
    this.rings = r.rings;
    this.atomRings = r.atomRings;
    this.bondRings = r.bondRings;
    this.aromAtoms = r.aromatic.atoms;
    this.aromBonds = r.aromatic.bonds;
    const seen = new Set();
    this.aromRings = [];
    for (const ring of r.aromatic.aromaticRings) {
      const k = [...ring].sort((a, b) => a - b).join(',');
      if (!seen.has(k)) { seen.add(k); this.aromRings.push(ring); }
    }
    this.sym = computeSymmetry(mol, this.aromBonds, this.atomRings);
    this.warnings = [];
    this._groups = new Map();
  }

  warn(msg) { if (!this.warnings.includes(msg)) this.warnings.push(msg); }
  atom(id) { return this.mol.getAtom(id); }
  el(id) { return this.mol.getAtom(id).el; }
  heavy(id) { return this.mol.neighbors(id).filter((n) => this.el(n) !== 'H'); }
  nH(id) { return this.mol.totalH(id); }
  isArom(id) { return this.aromAtoms.has(id); }
  bond(a, b) { return this.mol.bondBetween(a, b); }
  isRingBond(a, b) { const bd = this.bond(a, b); return !!bd && this.bondRings.has(bd.id); }

  /** Non-aromatic multiple-bond partners of `id` with the given order. */
  multi(id, order, except = null) {
    return this.mol.bondsOf(id)
      .filter((b) => b.order === order && !this.aromBonds.has(b.id))
      .map((b) => this.mol.otherAtom(b, id))
      .filter((n) => n !== except);
  }

  hyb(id) {
    if (this.isArom(id)) return 'ar';
    let d = 0, t = 0;
    for (const b of this.mol.bondsOf(id)) { if (b.order === 3) t++; else if (b.order === 2) d++; }
    if (t || d >= 2) return 'sp';
    return d ? 'sp2' : 'sp3';
  }

  isSp3C(id) { return this.el(id) === 'C' && this.hyb(id) === 'sp3'; }

  isCarbonylC(id) {
    return this.el(id) === 'C' && !this.isArom(id)
      && this.multi(id, 2).some((n) => this.el(n) === 'O' || this.el(n) === 'S');
  }

  isSulfonyl(id) {
    return (this.el(id) === 'S' || this.el(id) === 'P') && this.multi(id, 2).some((n) => this.el(n) === 'O');
  }

  sameRing(a, b) {
    const ra = this.atomRings.get(a) || [];
    const rb = this.atomRings.get(b) || [];
    return ra.some((x) => rb.includes(x));
  }

  smallestRing(id) {
    const idx = this.atomRings.get(id) || [];
    return idx.length ? Math.min(...idx.map((i) => this.rings[i].length)) : 0;
  }

  aromRingOf(id) {
    const c = this.aromRings.filter((r) => r.includes(id));
    if (!c.length) return null;
    return c.find((r) => r.length === 6) ?? c.sort((a, b) => a.length - b.length)[0];
  }

  /** π-donor heteroatom of a five-ring (pyrrole N, furan O, thiophene S). */
  isPiDonor(id) {
    const e = this.el(id);
    if (e === 'O' || e === 'S' || e === 'Se') return true;
    if (e === 'N') return this.mol.degree(id) + this.mol.implicitH(id) === 3 && this.atom(id).charge === 0;
    return false;
  }

  /** Substituent group key for atom `s` as seen from atom `from`. */
  group(from, s) {
    const key = `${from}>${s}`;
    if (!this._groups.has(key)) this._groups.set(key, this._classify(from, s));
    return this._groups.get(key);
  }

  _classify(from, s) {
    const e = this.el(s);
    if (e === 'H') return 'H';
    if (HALOGENS.has(e)) return e;
    const bd = this.bond(from, s);
    if (bd && bd.order === 2 && !this.aromBonds.has(bd.id)) {
      if (e === 'O' || e === 'S') return '=O';
      if (e === 'N') return 'N=X';
      if (e === 'C') return 'vinyl';
    }
    if (this.isArom(s)) {
      if (this.isArom(from) && this.isRingBond(from, s)) {
        const ring = this.aromRings.find((r) => r.includes(from) && r.includes(s));
        if (ring && ring.length === 5 && ring.some((r) => ['N', 'O'].includes(this.el(r)) && this.isPiDonor(r))) {
          return e === 'C' ? 'fused5' : 'fused5X';
        }
        return 'fused';
      }
      return e === 'N' ? 'Narom' : 'aryl';
    }
    switch (e) {
      case 'C': return this._carbonGroup(from, s);
      case 'N': return this._nitrogenGroup(from, s);
      case 'O': return this._oxygenGroup(from, s);
      case 'S': {
        const nO = this.multi(s, 2).filter((n) => this.el(n) === 'O').length;
        if (nO >= 2) return 'SO2R';
        if (nO === 1) return 'SOR';
        return this.nH(s) > 0 ? 'SH' : 'SR';
      }
      case 'Si': return 'SiR3';
      case 'P': return 'P';
      case 'B': return 'B';
      default: return 'other';
    }
  }

  _carbonGroup(from, s) {
    if (this.multi(s, 3, from).length) {
      return this.el(this.multi(s, 3, from)[0]) === 'N' ? 'CN' : 'alkynyl';
    }
    const dbl = this.multi(s, 2, from).map((n) => this.el(n));
    if (dbl.includes('O')) return this._carbonylType(s, from);
    if (dbl.includes('S')) return 'COR';
    if (dbl.includes('N')) return 'C=N';
    if (dbl.includes('C')) return 'vinyl';
    const hv = this.heavy(s).filter((n) => n !== from);
    const hal = hv.filter((n) => HALOGENS.has(this.el(n))).length;
    if (hal >= 2) return 'CX3';
    if (hv.some((n) => this.el(n) === 'O')) return 'alkylO';
    if (hv.some((n) => this.el(n) === 'N')) return 'alkylN';
    if (hal) return 'alkylHal';
    return ['CH3', 'alkyl', 'alkylCH', 'alkylC'][Math.min(hv.length, 3)];
  }

  _carbonylType(c, from) {
    if (this.nH(c) > 0) return 'CHO';
    const t = this.heavy(c).find((n) => n !== from && !(this.el(n) === 'O' && this.bond(c, n).order === 2));
    if (t === undefined) return 'COR';
    switch (this.el(t)) {
      case 'O':
        if (this.nH(t) > 0) return 'COOH';
        return this.atom(t).charge < 0 ? 'COO-' : 'COOR';
      case 'N': return 'CONR2';
      case 'Cl': case 'Br': case 'F': case 'I': return 'COCl';
      case 'S': return 'COOR';
      default: return 'COR';
    }
  }

  _nitrogenGroup(from, s) {
    const hv = this.heavy(s).filter((n) => n !== from);
    if (hv.filter((n) => this.el(n) === 'O').length >= 2) return 'NO2';
    if (this.multi(s, 2).length || this.multi(s, 3).length) return 'N=X';
    if (this.atom(s).charge > 0) return 'NR3+';
    if (hv.some((t) => this.isCarbonylC(t) || this.isSulfonyl(t))) return 'Namide';
    if (hv.some((t) => this.isArom(t))) return 'NAr';
    const h = this.nH(s);
    return h >= 2 ? 'NH2' : h === 1 ? 'NHR' : 'NR2';
  }

  _oxygenGroup(from, s) {
    if (this.atom(s).charge < 0) return 'O-';
    if (this.nH(s) > 0) return 'OH';
    const t = this.heavy(s).find((n) => n !== from);
    if (t === undefined) return 'OR';
    if (this.isCarbonylC(t)) {
      return this.heavy(t).some((n) => n !== s && this.isArom(n)) ? 'OCOAr' : 'OCOR';
    }
    if (this.isSulfonyl(t)) return 'OSO2';
    if (this.isArom(t)) return 'OAr';
    if (this.el(t) === 'C' && this.multi(t, 2).some((n) => this.el(n) === 'C')) return 'Ovinyl';
    return 'OR';
  }

  /**
   * Aromatic substituent walk: calls fn(groupKey, ringDistance, ringAtom) for every
   * non-ring heavy neighbour of every atom of `ring` (distance measured from `c`).
   */
  forRingSubstituents(ring, c, fn) {
    const n = ring.length;
    const ic = ring.indexOf(c);
    const inRing = new Set(ring);
    ring.forEach((r, i) => {
      const d = Math.min(Math.abs(i - ic), n - Math.abs(i - ic));
      for (const s of this.heavy(r)) if (!inRing.has(s)) fn(this.group(r, s), d, r);
    });
  }

  ringDist(ring, a, b) {
    const i = ring.indexOf(a), j = ring.indexOf(b), n = ring.length;
    return Math.min(Math.abs(i - j), n - Math.abs(i - j));
  }
}

// ===========================================================================
// 1H prediction
// ===========================================================================

/** Hydrogen objects (one per H) with alkene geometry resolved. */
function buildHydrogens(ctx) {
  const { mol } = ctx;
  const hs = [];
  const byParent = new Map();
  for (const a of mol.atoms.values()) {
    if (a.el === 'H' && mol.neighbors(a.id).some((n) => ctx.el(n) !== 'H')) continue;
    if (a.el === 'H') continue; // isolated H / H2: ignored
    const list = [];
    const explicit = mol.neighbors(a.id).filter((n) => ctx.el(n) === 'H');
    const total = mol.totalH(a.id);
    for (let k = 0; k < total; k++) {
      const h = { uid: `${a.id}:${k}`, parent: a.id, explicit: explicit[k] ?? null, side: 0, db: null };
      list.push(h);
      hs.push(h);
    }
    byParent.set(a.id, list);
  }
  ctx.hByParent = byParent;
  resolveAlkeneGeometry(ctx);
  // Keys: parent symmetry class, plus cis-partner label for inequivalent =CH2 protons.
  for (const h of hs) {
    let key = `${ctx.sym.get(h.parent)}`;
    if (h.db && byParent.get(h.parent).length === 2) {
      const cis = h.db.other.filter((sl) => sl.side === h.side)
        .map((sl) => (sl.h ? 'H' : `A${ctx.sym.get(sl.id)}`)).sort().join(',');
      key += `|cis:${cis}`;
    }
    h.key = key;
  }
  return hs;
}

/**
 * For every isolated C=C, assign each substituent / H a side (+1/−1) of the
 * double-bond axis, from 2D coordinates where available, else by ring
 * constraints (cis in small rings) or assuming the E isomer.
 */
function resolveAlkeneGeometry(ctx) {
  const { mol } = ctx;
  for (const b of mol.bonds.values()) {
    if (b.order !== 2 || ctx.aromBonds.has(b.id)) continue;
    const c1 = b.a1, c2 = b.a2;
    if (ctx.el(c1) !== 'C' || ctx.el(c2) !== 'C' || ctx.hyb(c1) !== 'sp2' || ctx.hyb(c2) !== 'sp2') continue;
    const p1 = ctx.atom(c1), p2 = ctx.atom(c2);
    const dx = p2.x - p1.x, dy = p2.y - p1.y;
    const len = Math.hypot(dx, dy);
    const geomSide = (id) => {
      if (id === null || id === undefined || len < 1e-6) return 0;
      const q = ctx.atom(id);
      const cr = dx * (q.y - p1.y) - dy * (q.x - p1.x);
      return Math.abs(cr) < 1e-6 * len ? 0 : Math.sign(cr);
    };
    const slots = (c, other) => [
      ...ctx.heavy(c).filter((n) => n !== other).map((id) => ({ id, side: geomSide(id) })),
      ...ctx.hByParent.get(c).map((h) => ({ h, id: h.explicit, side: geomSide(h.explicit) })),
    ];
    const A = slots(c1, c2), B = slots(c2, c1);
    const complete = (E) => {
      if (E.length !== 2) return;
      const [x, y] = E;
      if (x.side && !y.side) y.side = -x.side;
      else if (y.side && !x.side) x.side = -y.side;
      else if (x.side && x.side === y.side) y.side = -x.side; // inconsistent drawing
    };
    const known = (E) => E.some((s) => s.side !== 0);
    const pick = (E) => E.find((s) => !s.h && (!known(E) || s.side)) ?? E.find((s) => !known(E) || s.side);
    const fix = (U, K) => {
      const u = U.find((s) => !s.h) ?? U[0];
      const k = pick(K);
      if (!u) return;
      if (!k) u.side = 1;
      else if (!u.h && !k.h) u.side = ctx.sameRing(u.id, k.id) ? k.side : -k.side; // ring => cis, else E
      else u.side = 1;
      complete(U);
    };
    complete(A); complete(B);
    if (!known(A) && !known(B) && A.length) { (A.find((s) => !s.h) ?? A[0]).side = 1; complete(A); }
    if (!known(A)) fix(A, B);
    else if (!known(B)) fix(B, A);
    for (const s of A) if (s.h) { s.h.side = s.side; s.h.db = { bond: b, partner: c2, own: A, other: B }; }
    for (const s of B) if (s.h) { s.h.side = s.side; s.h.db = { bond: b, partner: c1, own: B, other: A }; }
  }
}

/** Shift + environment + exchangeable flag for one hydrogen. */
function hShift(ctx, h) {
  const c = h.parent;
  const e = ctx.el(c);
  if (e !== 'C') return exchangeableShift(ctx, c);
  if (ctx.isArom(c)) return aromaticHShift(ctx, c);
  const hy = ctx.hyb(c);
  if (hy === 'sp3') return sp3HShift(ctx, c);
  if (hy === 'sp') {
    const t = ctx.multi(c, 3)[0];
    if (t !== undefined && ctx.el(t) === 'C') {
      const g = ctx.heavy(t).find((n) => n !== c);
      const inc = g === undefined ? 0 : lookup(ALKYNE_H, ctx.group(t, g)) ?? 0;
      return { shift: ALKYNE_H_BASE + inc, env: '≡C–H (terminal alkyne)' };
    }
    return { shift: 4.7, env: 'allenic / cumulated CH' };
  }
  // sp2
  const dbl = ctx.multi(c, 2)[0];
  const de = ctx.el(dbl);
  if (de === 'O' || de === 'S') {
    const t = ctx.heavy(c).find((n) => n !== dbl);
    if (t === undefined) return { shift: 9.6, env: 'formaldehyde' };
    const te = ctx.el(t);
    if (te === 'O') return { shift: 8.05, env: 'HC(=O)O (formate)' };
    if (te === 'N') return { shift: 8.05, env: 'HC(=O)N (formamide)' };
    const g = ctx.group(c, t);
    if (g === 'aryl' || g === 'fused') return { shift: 9.95, env: 'CHO (aromatic aldehyde)' };
    if (g === 'vinyl') return { shift: 9.5, env: 'CHO (conjugated aldehyde)' };
    return { shift: 9.75, env: 'CHO (aldehyde)' };
  }
  if (de === 'N') {
    const oxime = ctx.heavy(dbl).some((n) => ctx.el(n) === 'O');
    return { shift: oxime ? 7.5 : 8.2, env: oxime ? 'CH=N–O (oxime)' : 'CH=N (imine)' };
  }
  if (de === 'C' && ctx.multi(dbl, 2).length >= 2) return { shift: 4.7, env: '=CH2 (allene)' };
  if (de === 'C' && h.db) return vinylHShift(ctx, h);
  return { shift: 5.3, env: 'vinyl CH' };
}

function sp3HShift(ctx, c) {
  const nh = Math.min(ctx.nH(c), 4);
  const base = H_SP3_BASE[nh];
  const ti = nh >= 3 ? 0 : nh === 2 ? 1 : 2;
  const alpha = [], beta = [], labels = [];
  for (const s of ctx.heavy(c)) {
    const g = ctx.group(c, s);
    if (ALKYL.has(g)) {
      if (!ctx.isSp3C(s)) continue;
      for (const t of ctx.heavy(s)) {
        if (t === c) continue;
        const gb = ctx.group(s, t);
        if (ALKYL.has(gb)) continue;
        const v = lookup(H_SP3_BETA, gb);
        if (v !== null) beta.push(v);
      }
    } else {
      const inc = lookup(H_SP3_ALPHA, g);
      if (inc) alpha.push(inc[ti]);
      labels.push(label(g));
    }
  }
  let shift = base + damped(alpha) + damped(beta);
  const ring = ctx.smallestRing(c);
  if (ring === 3) shift -= 1.0;
  else if (ring === 4) shift += 0.5;
  const name = ['CH', 'CH2', 'CH3', 'CH4'][nh - 1];
  const env = labels.length ? `${name} next to ${[...new Set(labels)].join(', ')}` : `${name} (alkyl)`;
  return { shift: clamp(shift, 0, 7.5), env };
}

function aromaticHShift(ctx, c) {
  const ring = ctx.aromRingOf(c);
  if (!ring) return { shift: 7.26, env: 'aromatic CH' };
  let shift = 7.26;
  let env = 'aromatic CH';
  const subs = [];
  if (ring.length === 5) {
    const donors = ring.filter((r) => r !== c && ctx.isPiDonor(r));
    if (donors.length) {
      const x = donors.sort((a, b) => ctx.ringDist(ring, c, a) - ctx.ringDist(ring, c, b))[0];
      const pos = ctx.ringDist(ring, c, x) === 1 ? 'a' : 'b';
      const tbl = FIVE_RING_H[ctx.el(x)] ?? FIVE_RING_H.S;
      shift = tbl[pos];
      env = `heteroaromatic CH (${pos === 'a' ? 'α' : 'β'} to ${ctx.el(x)})`;
      for (const r of ring) {
        if (r === c || r === x || ctx.el(r) !== 'N' || ctx.isPiDonor(r)) continue;
        shift += AZOLE_N_H[Math.min(ctx.ringDist(ring, c, r), 2) - 1];
      }
    } else {
      shift = 5.7;
      env = 'five-ring aromatic CH';
    }
    ctx.forRingSubstituents(ring, c, (g, d) => {
      const inc = lookup(AROMATIC_H, g);
      if (inc && d >= 1) {
        shift += inc[Math.min(d, 2) - 1];
        if (!g.startsWith('fused')) subs.push(`${POS_NAME[Math.min(d, 2)]} to ${label(g)}`);
      }
    });
  } else {
    let nAlpha = 0;
    const pos = [];
    for (const r of ring) {
      if (r === c || ctx.el(r) === 'C') continue;
      if (ctx.el(r) !== 'N') { ctx.warn(`Unusual aromatic ring atom ${ctx.el(r)}; estimate is rough`); continue; }
      const d = Math.min(ctx.ringDist(ring, c, r), 3);
      if (d === 1) shift += PYRIDINE_N_H[0] * (nAlpha++ ? ALPHA_N_DAMP_H : 1);
      else shift += PYRIDINE_N_H[d - 1];
      pos.push(['α', 'β', 'γ'][d - 1]);
    }
    if (pos.length) env = `heteroaromatic CH (${pos.join(', ')} to N)`;
    ctx.forRingSubstituents(ring, c, (g, d) => {
      const inc = lookup(AROMATIC_H, g);
      if (inc && d >= 1) { shift += inc[d - 1]; if (!g.startsWith('fused')) subs.push(`${POS_NAME[d]} to ${label(g)}`); }
    });
  }
  if (subs.length) env += ` (${subs.join(', ')})`;
  return { shift: clamp(shift, 5.5, 9.8), env };
}

function vinylHShift(ctx, h) {
  const c = h.parent;
  const { partner, own, other } = h.db;
  let shift = 5.25;
  const parts = [];
  for (const sl of own) {
    if (sl.h) continue;
    const g = ctx.group(c, sl.id);
    const z = lookup(VINYL_H, g);
    if (z) shift += z[0];
    parts.push(`gem ${label(g)}`);
  }
  for (const sl of other) {
    if (sl.h) continue;
    const g = ctx.group(partner, sl.id);
    const z = lookup(VINYL_H, g);
    const cis = sl.side === h.side;
    if (z) shift += cis ? z[1] : z[2];
    parts.push(`${cis ? 'cis' : 'trans'} ${label(g)}`);
  }
  const nh = ctx.hByParent.get(c).length;
  const env = `${nh === 2 ? '=CH2' : '=CH'} (vinyl${parts.length ? `; ${parts.join(', ')}` : ''})`;
  return { shift: clamp(shift, 3.8, 8.0), env };
}

function exchangeableShift(ctx, x) {
  const e = ctx.el(x);
  const hv = ctx.heavy(x);
  const res = (shift, env, exchangeable = true) => ({ shift, env, exchangeable });
  if (e === 'O') {
    const t = hv[0];
    if (t === undefined) return res(1.56, 'H2O');
    if (ctx.isCarbonylC(t)) {
      const conj = ctx.heavy(t).some((n) => n !== x && CONJ.has(ctx.group(t, n)));
      return res(conj ? 12.0 : 11.5, 'COOH (carboxylic acid)');
    }
    if (ctx.isArom(t)) {
      const ring = ctx.aromRingOf(t);
      let hbond = false;
      if (ring) {
        for (const r of ring) {
          if (ctx.ringDist(ring, t, r) !== 1) continue;
          for (const s of ctx.heavy(r)) {
            if (!ring.includes(s) && ['CHO', 'COR', 'COOR', 'COOH', 'CONR2', 'NO2'].includes(ctx.group(r, s))) hbond = true;
          }
        }
      }
      return res(hbond ? 10.8 : 5.0, hbond ? 'phenol OH (H-bonded)' : 'phenol OH');
    }
    if (ctx.isSulfonyl(t)) return res(10.0, 'acidic OH (S/P–OH)');
    if (ctx.el(t) === 'N') return res(8.5, 'N–OH (oxime / hydroxylamine)');
    if (ctx.el(t) === 'C' && ctx.hyb(t) === 'sp2') return res(5.5, 'enol OH');
    if (ctx.el(t) !== 'C') return res(4.5, 'OH');
    return res(2.0, 'OH (alcohol)');
  }
  if (e === 'N') {
    if (ctx.atom(x).charge > 0) return res(7.5, 'N+–H (ammonium)');
    if (ctx.isArom(x)) return res(8.0, 'NH (pyrrole-type)');
    const carbonyls = hv.filter((t) => ctx.isCarbonylC(t));
    if (carbonyls.length >= 2) return res(8.5, 'NH (imide)');
    if (carbonyls.length === 1) {
      const t = carbonyls[0];
      const het = ctx.heavy(t).filter((n) => n !== x && ctx.bond(t, n).order === 1).map((n) => ctx.el(n));
      if (hv.some((n) => ctx.isArom(n))) return res(7.6, 'NH (anilide)');
      if (het.includes('O')) return res(4.9, 'NH (carbamate)');
      if (het.includes('N')) return res(5.0, 'NH (urea)');
      if (ctx.nH(x) >= 2) return res(6.0, 'NH2 (primary amide)');
      return res(5.9, 'NH (amide)');
    }
    if (hv.some((t) => ctx.isSulfonyl(t))) return res(4.8, 'NH (sulfonamide)');
    if (ctx.multi(x, 2).length) return res(8.0, '=NH (imine)');
    if (hv.some((t) => ctx.isArom(t))) return res(3.6, ctx.nH(x) >= 2 ? 'NH2 (aniline)' : 'NH (N-aryl amine)');
    if (hv.some((t) => ['N', 'O'].includes(ctx.el(t)))) return res(3.5, 'NH (hydrazine / hydroxylamine)');
    return res(1.3, ctx.nH(x) >= 2 ? 'NH2 (amine)' : 'NH (amine)');
  }
  if (e === 'S') return res(hv.some((t) => ctx.isArom(t)) ? 3.4 : 1.4, 'SH (thiol)');
  if (e === 'Si') return res(3.9, 'Si–H', false);
  if (e === 'B') return res(1.0, 'B–H', false);
  if (e === 'P') return res(3.5, 'P–H', false);
  ctx.warn(`No 1H model for H on ${e}; shift is a placeholder`);
  return res(2.0, `${e}–H`, false);
}

/** Vicinal (and =CH2 geminal) coupling constant between H objects h (on c) and p (on n). */
function couplingJ(ctx, h, c, p, n) {
  if (c === n) return J.gem;
  const bd = ctx.bond(c, n);
  if (ctx.isArom(c) && ctx.isArom(n) && ctx.aromBonds.has(bd.id)) {
    const ring = ctx.aromRingOf(c);
    if (ring && ring.length === 5) return J.fiveRing;
    const nextToN = [c, n].some((a) => ctx.heavy(a).some((x) => ctx.el(x) === 'N' && ctx.isArom(x)));
    return nextToN ? J.orthoNextToN : J.ortho;
  }
  if (bd.order === 2) return h.side && h.side === p.side ? J.cis : J.trans;
  const hc = ctx.hyb(c), hn = ctx.hyb(n);
  if (ctx.isCarbonylC(c) || ctx.isCarbonylC(n)) {
    const otherSp3 = ctx.isCarbonylC(c) ? hn === 'sp3' : hc === 'sp3';
    return otherSp3 ? J.aldehydeSp3 : J.aldehydeSp2;
  }
  if (hc !== 'sp3' && hn !== 'sp3') return J.diene;
  if (hc !== 'sp3' || hn !== 'sp3') return J.allylic;
  return J.alkyl;
}

function binomial(n) {
  const row = [1];
  for (let k = 1; k <= n; k++) row.push((row[k - 1] * (n - k + 1)) / k);
  return row;
}

/** First-order line list (Pascal's triangle) for one signal. */
function buildPeaks(shift, nH, groups, freq) {
  let lines = [{ off: 0, w: 1 }];
  for (const { J: j, n } of groups) {
    const coeffs = binomial(n);
    const tot = 2 ** n;
    const next = [];
    for (const l of lines) for (let k = 0; k <= n; k++) next.push({ off: l.off + (k - n / 2) * j, w: (l.w * coeffs[k]) / tot });
    next.sort((a, b) => a.off - b.off);
    lines = [];
    for (const l of next) {
      const last = lines[lines.length - 1];
      if (last && Math.abs(last.off - l.off) < 0.01) last.w += l.w;
      else lines.push({ ...l });
    }
  }
  return lines.map((l) => ({ shift: round(shift + l.off / freq, 5), intensity: l.w * nH }));
}

/**
 * Predict the 1H NMR spectrum.
 * @param {import('./molecule.js').Molecule} mol
 * @param {{frequency?: number}} [opts] spectrometer frequency in MHz (affects peaks / second-order flags)
 * @returns {{signals: Array<{shift:number, atomIds:number[], nH:number, multiplicity:string, J:number[],
 *   environment:string, exchangeable:boolean, peaks:Array<{shift:number,intensity:number}>}>, warnings:string[]}}
 */
export function predictH1(mol, { frequency = 400 } = {}) {
  const ctx = new Ctx(mol);
  const hs = buildHydrogens(ctx);
  const info = new Map(); // key -> {shift, env, exchangeable}
  const groups = new Map(); // key -> H objects
  for (const h of hs) {
    if (!groups.has(h.key)) {
      const r = hShift(ctx, h);
      info.set(h.key, { shift: round(r.shift, 2), env: r.env, exchangeable: !!r.exchangeable });
      groups.set(h.key, []);
    }
    groups.get(h.key).push(h);
  }
  const radicals = mol.atomList().some((a) => a.radical);
  if (radicals) ctx.warn('Radicals present: predictions are not meaningful for paramagnetic species');

  const signals = [];
  for (const [key, list] of groups) {
    const { shift, env, exchangeable } = info.get(key);
    const rep = list[0];
    const nH = list.length;
    const atomIds = [...new Set(list.map((h) => h.parent))].sort((a, b) => a - b);
    let multiplicity = 's', Js = [], jGroups = [];
    if (exchangeable) {
      multiplicity = 'br s';
    } else if (ctx.el(rep.parent) === 'C') {
      // Collect coupling partners of the representative proton.
      const partners = new Map(); // key -> {n, J}
      const add = (p, j) => {
        if (p.key === key || info.get(p.key)?.exchangeable) return;
        const cur = partners.get(p.key) ?? { n: 0, J: 0 };
        cur.J = (cur.J * cur.n + j) / (cur.n + 1);
        cur.n++;
        partners.set(p.key, cur);
      };
      for (const p of ctx.hByParent.get(rep.parent)) if (p !== rep) add(p, J.gem);
      for (const nb of ctx.heavy(rep.parent)) {
        if (ctx.el(nb) !== 'C') continue;
        for (const p of ctx.hByParent.get(nb)) add(p, couplingJ(ctx, rep, rep.parent, p, nb));
      }
      let secondOrder = false;
      for (const [pk, { J: j }] of partners) {
        if ((Math.abs(info.get(pk).shift - shift) * frequency) / j < STRONG_COUPLING_RATIO) secondOrder = true;
      }
      // Merge partner sets with (nearly) equal J (n+1 rule).
      const sorted = [...partners.values()].sort((a, b) => b.J - a.J);
      for (const p of sorted) {
        const last = jGroups[jGroups.length - 1];
        if (last && Math.abs(last.J - p.J) <= J_MERGE_TOL) {
          last.J = (last.J * last.n + p.J * p.n) / (last.n + p.n);
          last.n += p.n;
        } else jGroups.push({ ...p });
      }
      if (jGroups.length === 0) multiplicity = 's';
      else if (jGroups.length === 1) multiplicity = MULT_NAMES[jGroups[0].n] ?? 'm';
      else if (jGroups.length === 2) {
        const [a, b] = jGroups.map((g) => g.n);
        multiplicity = a === 1 && b === 1 ? 'dd' : a === 1 && b === 2 ? 'dt' : a === 2 && b === 1 ? 'td' : 'm';
      } else multiplicity = 'm';
      if (secondOrder || partners.size > 3) multiplicity = jGroups.length ? 'm' : multiplicity;
      if (multiplicity !== 's' && multiplicity !== 'm') Js = jGroups.map((g) => round(g.J, 1));
    }
    signals.push({
      shift, atomIds, nH, multiplicity, J: Js, environment: env, exchangeable,
      peaks: buildPeaks(shift, nH, exchangeable ? [] : jGroups, frequency),
    });
  }
  signals.sort((a, b) => b.shift - a.shift || b.nH - a.nH);
  return { signals, warnings: ctx.warnings };
}

// ===========================================================================
// 13C prediction
// ===========================================================================

function c13Shift(ctx, c) {
  if (ctx.atom(c).charge !== 0 || ctx.atom(c).radical) ctx.warn('Charged or radical carbon: 13C estimate is rough');
  if (ctx.isArom(c)) return c13Aromatic(ctx, c);
  const hy = ctx.hyb(c);
  if (hy === 'sp3') return c13Sp3(ctx, c);
  const tri = ctx.multi(c, 3);
  if (tri.length) {
    const t = tri[0];
    if (ctx.el(t) === 'N') {
      const g = ctx.heavy(c).find((n) => n !== t);
      const conj = g !== undefined && CONJ.has(ctx.group(c, g));
      return { shift: conj ? C13_CO.nitrileConj : C13_CO.nitrile, env: 'C≡N (nitrile)' };
    }
    return c13Alkyne(ctx, c, t);
  }
  const dbl = ctx.multi(c, 2);
  if (dbl.length >= 2) {
    const els = dbl.map((n) => ctx.el(n)).sort().join('');
    if (els === 'CC') return { shift: C13_CO.alleneCentral, env: '=C= (allene central)' };
    if (els === 'NO' || els === 'NS') return { shift: C13_CO.isocyanate, env: 'N=C=O/S' };
    if (els === 'OO') return { shift: C13_CO.co2, env: 'O=C=O' };
    if (els === 'CO') return { shift: C13_CO.ketene, env: 'C=C=O (ketene)' };
    ctx.warn('Unusual cumulated carbon');
    return { shift: 150, env: 'cumulated C' };
  }
  const d = dbl[0];
  const de = ctx.el(d);
  if (de === 'O' || de === 'S') return c13Carbonyl(ctx, c, d);
  if (de === 'N') return c13Imine(ctx, c, d);
  if (de === 'C') {
    if (ctx.multi(d, 2).length >= 2) return { shift: C13_CO.alleneTerminal, env: '=C (allene terminal)' };
    return c13Alkene(ctx, c, d);
  }
  ctx.warn(`No 13C model for C=${de}`);
  return { shift: 150, env: `C=${de}` };
}

function c13Sp3(ctx, c) {
  let shift = C13_SP3.base;
  const n = [0, 0, 0, 0, 0];
  const alphaInc = [], betaInc = [];
  const labels = [];
  const dist = new Map([[c, 0]]);
  const queue = [c];
  while (queue.length) {
    const cur = queue.shift();
    const d = dist.get(cur);
    for (const nb of ctx.heavy(cur)) {
      if (dist.has(nb)) continue;
      if (ctx.isSp3C(nb)) {
        if (d + 1 <= 4) {
          dist.set(nb, d + 1);
          n[d + 1]++;
          if (d + 1 < 4) queue.push(nb);
        }
      } else {
        dist.set(nb, d + 1);
        if (d + 1 > 3) continue;
        const g = ctx.group(cur, nb);
        const inc = lookup(C13_SP3_INC, g);
        if (!inc) continue;
        if (d === 0) {
          alphaInc.push(inc[0]);
          labels.push(label(g));
          const extra = C13_HALOGEN_BRANCH[g];
          if (extra) {
            const nc = ctx.heavy(c).filter((x) => ctx.el(x) === 'C').length;
            shift += extra * Math.max(0, nc - 1);
          }
        } else if (d === 1) betaInc.push(inc[1]);
        else shift += inc[2];
      }
    }
  }
  shift += C13_SP3.alpha * n[1] + C13_SP3.beta * n[2] + C13_SP3.gamma * n[3] + C13_SP3.delta * n[4];
  shift += damped(alphaInc) + damped(betaInc);
  // Grant–Paul branching corrections.
  const deg = (id) => ctx.heavy(id).filter((x) => ctx.el(x) === 'C').length;
  const di = deg(c);
  for (const nb of ctx.heavy(c)) {
    if (!ctx.isSp3C(nb)) continue;
    shift += C13_BRANCH[di]?.[deg(nb)] ?? 0;
  }
  shift += C13_RING[ctx.smallestRing(c)] ?? 0;
  const env = labels.length ? `sp3 C next to ${[...new Set(labels)].join(', ')}` : 'sp3 C (alkyl)';
  return { shift: clamp(shift, 0, 110), env };
}

function c13Aromatic(ctx, c) {
  const ring = ctx.aromRingOf(c);
  if (!ring) return { shift: 128.5, env: 'aromatic C' };
  let shift = 128.5;
  let env = 'aromatic C';
  if (ring.length === 5) {
    const donors = ring.filter((r) => r !== c && ctx.isPiDonor(r));
    if (donors.length) {
      const x = donors.sort((a, b) => ctx.ringDist(ring, c, a) - ctx.ringDist(ring, c, b))[0];
      const pos = ctx.ringDist(ring, c, x) === 1 ? 'a' : 'b';
      shift = (FIVE_RING_C[ctx.el(x)] ?? FIVE_RING_C.S)[pos];
      env = `heteroaromatic C (${pos === 'a' ? 'α' : 'β'} to ${ctx.el(x)})`;
      let first = true;
      for (const r of ring) {
        if (r === c || r === x || ctx.el(r) !== 'N' || ctx.isPiDonor(r)) continue;
        const d = Math.min(ctx.ringDist(ring, c, r), 2);
        shift += AZOLE_N_C[d - 1] * (d === 1 && !first ? 0.6 : 1);
        if (d === 1) first = false;
      }
    } else shift = 103;
  } else {
    let nAlpha = 0;
    const pos = [];
    for (const r of ring) {
      if (r === c || ctx.el(r) === 'C') continue;
      if (ctx.el(r) !== 'N') { ctx.warn(`Unusual aromatic ring atom ${ctx.el(r)}; estimate is rough`); continue; }
      const d = Math.min(ctx.ringDist(ring, c, r), 3);
      if (d === 1) shift += PYRIDINE_N_C[0] * (nAlpha++ ? ALPHA_N_DAMP_C : 1);
      else shift += PYRIDINE_N_C[d - 1];
      pos.push(['α', 'β', 'γ'][d - 1]);
    }
    if (pos.length) env = `heteroaromatic C (${pos.join(', ')} to N)`;
  }
  const maxD = ring.length === 5 ? 2 : 3;
  let ipso = null;
  ctx.forRingSubstituents(ring, c, (g, d) => {
    const inc = lookup(AROMATIC_C, g);
    if (inc) shift += inc[Math.min(d, maxD)];
    const fusion = g.startsWith('fused');
    if (d === 0 && !fusion) ipso = g;
    if (d === 0 && fusion && !ipso) ipso = 'ring fusion';
  });
  if (ipso) env += ` (ipso, ${label(ipso)})`;
  return { shift: clamp(shift, 95, 170), env };
}

/** Plain alkyl chain counts (β, γ carbons) for a substituent carbon g attached to c. */
function chainCounts(ctx, c, g) {
  const beta = ctx.heavy(g).filter((x) => x !== c && ctx.el(x) === 'C');
  let gamma = 0;
  for (const b of beta) gamma += ctx.heavy(b).filter((x) => x !== g && ctx.el(x) === 'C').length;
  return { nb: beta.length, ng: gamma };
}

function c13Alkene(ctx, c, c2) {
  let shift = ALKENE_C_BASE;
  const K = ALKENE_C_ALKYL;
  for (const g of ctx.heavy(c)) {
    if (g === c2) continue;
    const key = ctx.group(c, g);
    if (PLAIN_ALKYL.has(key)) { const { nb, ng } = chainCounts(ctx, c, g); shift += K.a + K.b * nb + K.g * ng; } else shift += lookup(ALKENE_C, key)?.[0] ?? 0;
  }
  for (const g of ctx.heavy(c2)) {
    if (g === c) continue;
    const key = ctx.group(c2, g);
    if (PLAIN_ALKYL.has(key)) { const { nb, ng } = chainCounts(ctx, c2, g); shift += K.a2 + K.b2 * nb + K.g2 * ng; } else shift += lookup(ALKENE_C, key)?.[1] ?? 0;
  }
  const nh = ctx.nH(c);
  return { shift: clamp(shift, 75, 175), env: nh === 2 ? '=CH2 (alkene)' : nh === 1 ? '=CH (alkene)' : '=C (alkene)' };
}

function c13Alkyne(ctx, c, c2) {
  let shift = ALKYNE_C_BASE;
  const K = ALKYNE_C_ALKYL;
  for (const g of ctx.heavy(c)) {
    if (g === c2) continue;
    const key = ctx.group(c, g);
    if (PLAIN_ALKYL.has(key)) shift += K.a + K.b * chainCounts(ctx, c, g).nb;
    else shift += lookup(ALKYNE_C, key)?.[0] ?? 0;
  }
  for (const g of ctx.heavy(c2)) {
    if (g === c) continue;
    const key = ctx.group(c2, g);
    if (PLAIN_ALKYL.has(key)) shift += K.a2 + K.b2 * chainCounts(ctx, c2, g).nb;
    else shift += lookup(ALKYNE_C, key)?.[1] ?? 0;
  }
  return { shift: clamp(shift, 55, 100), env: ctx.nH(c) ? '≡CH (alkyne)' : '≡C (alkyne)' };
}

function c13Carbonyl(ctx, c, o) {
  const T = C13_CO;
  const others = ctx.heavy(c).filter((n) => n !== o);
  const els = others.map((n) => ctx.el(n));
  const nO = els.filter((e) => e === 'O').length;
  const nN = els.filter((e) => e === 'N').length;
  const carbons = others.filter((n) => ctx.el(n) === 'C');
  const conj = carbons.some((n) => CONJ.has(ctx.group(c, n)));
  let nBeta = 0;
  for (const t of carbons) if (ctx.isSp3C(t)) nBeta += ctx.heavy(t).filter((x) => x !== c && ctx.el(x) === 'C').length;
  const ring = ctx.smallestRing(c);
  const ringCorr = ring === 5 ? T.ring5 : ring === 4 ? T.ring4 : 0;
  const beta3 = T.betaPerC * Math.min(nBeta, 3);

  if (ctx.el(o) === 'S') {
    if (nN >= 2) return { shift: T.thiourea, env: 'C=S (thiourea)' };
    if (nN === 1) return { shift: T.thioamide, env: 'C=S (thioamide)' };
    ctx.warn('Thiocarbonyl estimate is rough');
    return { shift: T.thione, env: 'C=S' };
  }
  if (nO + nN >= 2) {
    if (nO >= 2) return { shift: T.carbonate, env: 'C=O (carbonate)' };
    if (nN >= 2) return { shift: T.urea, env: 'C=O (urea)' };
    return { shift: T.carbamate, env: 'C=O (carbamate)' };
  }
  if (ctx.nH(c) > 0) {
    if (nO) return { shift: T.formate, env: 'HC=O (formate)' };
    if (nN) return { shift: T.formamide, env: 'HC=O (formamide)' };
    if (!carbons.length) return { shift: T.formaldehyde, env: 'formaldehyde' };
    const g = ctx.group(c, carbons[0]);
    if (g === 'aryl' || g === 'fused') return { shift: T.aldehydeConj, env: 'CHO (aromatic aldehyde)' };
    if (g === 'vinyl') return { shift: T.aldehydeVinyl, env: 'CHO (conjugated aldehyde)' };
    return { shift: T.aldehyde + T.aldehydeBeta * Math.min(nBeta, 3), env: 'CHO (aldehyde)' };
  }
  if (nO === 1) {
    const ox = others[els.indexOf('O')];
    if (ctx.nH(ox) > 0) return { shift: conj ? T.acidConj : T.acid + beta3, env: 'C=O (carboxylic acid)' };
    if (ctx.atom(ox).charge < 0) return { shift: T.carboxylate, env: 'CO2− (carboxylate)' };
    if (ctx.heavy(ox).some((n) => n !== c && ctx.isCarbonylC(n))) {
      return { shift: conj ? T.anhydrideConj : T.anhydride, env: 'C=O (anhydride)' };
    }
    return { shift: (conj ? T.esterConj : T.ester + beta3) + ringCorr, env: ring ? 'C=O (lactone)' : 'C=O (ester)' };
  }
  if (nN === 1) return { shift: (conj ? T.amideConj : T.amide + beta3) + ringCorr, env: ring ? 'C=O (lactam)' : 'C=O (amide)' };
  if (els.some((e) => HALOGENS.has(e))) return { shift: conj ? T.acidHalideConj : T.acidHalide, env: 'C=O (acid halide)' };
  if (els.includes('S')) return { shift: conj ? T.thioesterConj : T.thioester, env: 'C=O (thioester)' };
  // Ketone
  let shift = T.ketone + T.ketoneBeta * Math.min(nBeta, 4) + ringCorr;
  let nConj = 0;
  for (const t of carbons) {
    const g = ctx.group(c, t);
    if (CONJ.has(g)) shift += nConj++ ? T.ketoneConj2 : T.ketoneConj;
    else if (g === 'alkynyl') shift += T.ketoneAlkynyl;
    else if (['COR', 'CHO', 'COOR', 'COOH'].includes(g)) shift += T.ketoneDione;
  }
  return { shift: clamp(shift, 175, 225), env: 'C=O (ketone)' };
}

function c13Imine(ctx, c, n) {
  const T = C13_CO;
  const others = ctx.heavy(c).filter((x) => x !== n);
  const singleN = others.filter((x) => ctx.el(x) === 'N').length;
  const singleO = others.some((x) => ctx.el(x) === 'O');
  const onN = ctx.heavy(n).filter((x) => x !== c).map((x) => ctx.el(x));
  const conj = others.some((x) => ctx.el(x) === 'C' && CONJ.has(ctx.group(c, x)));
  if (singleN >= 2) return { shift: T.guanidine, env: 'C=N (guanidine)' };
  if (singleN === 1) return { shift: T.amidine, env: 'C=N (amidine)' };
  if (singleO) return { shift: T.imidate, env: 'C=N (imidate)' };
  const h = ctx.nH(c) > 0;
  if (onN.includes('O')) return { shift: h ? T.aldoxime : T.ketoxime, env: 'C=N–O (oxime)' };
  if (onN.includes('N')) return { shift: T.hydrazone, env: 'C=N–N (hydrazone)' };
  if (h) return { shift: conj ? T.aldimineConj : T.aldimine, env: 'CH=N (imine)' };
  return { shift: conj ? T.ketimineConj : T.ketimine, env: 'C=N (imine)' };
}

/**
 * Predict the proton-decoupled 13C NMR spectrum.
 * @param {import('./molecule.js').Molecule} mol
 * @returns {{signals: Array<{shift:number, atomIds:number[], environment:string,
 *   multiplicity:'C'|'CH'|'CH2'|'CH3', count:number}>, warnings:string[]}}
 */
export function predictC13(mol) {
  const ctx = new Ctx(mol);
  const byClass = new Map();
  for (const a of mol.atoms.values()) {
    if (a.el !== 'C') continue;
    const k = ctx.sym.get(a.id);
    if (!byClass.has(k)) byClass.set(k, []);
    byClass.get(k).push(a.id);
  }
  const signals = [];
  for (const ids of byClass.values()) {
    const r = c13Shift(ctx, ids[0]);
    const nh = ctx.nH(ids[0]);
    signals.push({
      shift: round(r.shift, 1),
      atomIds: ids.sort((a, b) => a - b),
      environment: r.env,
      multiplicity: ['C', 'CH', 'CH2', 'CH3'][Math.min(nh, 3)],
      count: ids.length,
    });
  }
  signals.sort((a, b) => b.shift - a.shift);
  return { signals, warnings: ctx.warnings };
}

// ===========================================================================
// Shared utilities
// ===========================================================================

/**
 * Topological symmetry classes (Morgan-style refinement over element, charge,
 * degree, H count, aromaticity and bond orders).
 * @returns {Map<number, number>} atomId -> class index (equal index = equivalent)
 */
export function symmetryClasses(mol) {
  const r = analyzeRings(mol);
  return computeSymmetry(mol, r.aromatic.bonds, r.atomRings);
}

/**
 * Lorentzian line-shape simulation.
 * @param {Array<{shift:number, intensity:number}>} peaks
 * @param {{min?:number, max?:number, points?:number, lineWidth?:number}} [opts] lineWidth = HWHM in ppm
 * @returns {{x: Float64Array, y: Float64Array}} x ascending (ppm), y normalised to max 1
 */
export function simulateSpectrum(peaks, { min, max, points = 4000, lineWidth = 0.004 } = {}) {
  const shifts = peaks.map((p) => p.shift);
  const lo = min ?? (shifts.length ? Math.min(...shifts) - 0.5 : 0);
  const hi = max ?? (shifts.length ? Math.max(...shifts) + 0.5 : 10);
  const x = new Float64Array(points);
  const y = new Float64Array(points);
  const step = points > 1 ? (hi - lo) / (points - 1) : 0;
  const w2 = lineWidth * lineWidth;
  for (let i = 0; i < points; i++) x[i] = lo + i * step;
  for (const { shift, intensity } of peaks) {
    for (let i = 0; i < points; i++) {
      const d = x[i] - shift;
      y[i] += (intensity * w2) / (d * d + w2);
    }
  }
  let m = 0;
  for (let i = 0; i < points; i++) if (y[i] > m) m = y[i];
  if (m > 0) for (let i = 0; i < points; i++) y[i] /= m;
  return { x, y };
}

/**
 * Journal-style 1H report, e.g.
 * "1H NMR (400 MHz, CDCl3) δ 7.26 (m, 5H), 3.72 (q, J = 7.1 Hz, 2H), 1.25 (t, J = 7.1 Hz, 3H)."
 * Adjacent overlapping multiplets ('m', < 0.2 ppm apart) are reported as a range.
 */
export function formatH1Report(result, { frequency = 400, solvent = 'CDCl3' } = {}) {
  const sig = [...result.signals].sort((a, b) => b.shift - a.shift);
  const items = [];
  for (let i = 0; i < sig.length; i++) {
    const s = sig[i];
    if (s.multiplicity === 'm') {
      let j = i, nH = s.nH;
      while (j + 1 < sig.length && sig[j + 1].multiplicity === 'm' && sig[j].shift - sig[j + 1].shift < 0.2) {
        j++;
        nH += sig[j].nH;
      }
      const range = j > i && s.shift !== sig[j].shift ? `${s.shift.toFixed(2)}–${sig[j].shift.toFixed(2)}` : s.shift.toFixed(2);
      items.push(`${range} (m, ${nH}H)`);
      i = j;
      continue;
    }
    const jtxt = s.J && s.J.length ? `J = ${s.J.map((v) => v.toFixed(1)).join(', ')} Hz, ` : '';
    items.push(`${s.shift.toFixed(2)} (${s.multiplicity}, ${jtxt}${s.nH}H)`);
  }
  return `1H NMR (${frequency} MHz, ${solvent}) δ ${items.length ? items.join(', ') : 'no signals'}.`;
}

/** Journal-style 13C report, e.g. "13C NMR (100 MHz, CDCl3) δ 137.9, 129.1, 21.5." */
export function formatC13Report(result, { frequency = 100, solvent = 'CDCl3' } = {}) {
  const items = [...result.signals].sort((a, b) => b.shift - a.shift).map((s) => s.shift.toFixed(1));
  return `13C NMR (${frequency} MHz, ${solvent}) δ ${items.length ? items.join(', ') : 'no signals'}.`;
}
