// Periodic table data used throughout Benzene.
// Masses are IUPAC standard atomic weights (abridged). Isotope data is provided
// for elements commonly encountered in organic chemistry (used for exact mass
// and isotope pattern calculations).

const RAW = `1 H Hydrogen 1.008
2 He Helium 4.0026
3 Li Lithium 6.94
4 Be Beryllium 9.0122
5 B Boron 10.81
6 C Carbon 12.011
7 N Nitrogen 14.007
8 O Oxygen 15.999
9 F Fluorine 18.998
10 Ne Neon 20.180
11 Na Sodium 22.990
12 Mg Magnesium 24.305
13 Al Aluminium 26.982
14 Si Silicon 28.085
15 P Phosphorus 30.974
16 S Sulfur 32.06
17 Cl Chlorine 35.45
18 Ar Argon 39.95
19 K Potassium 39.098
20 Ca Calcium 40.078
21 Sc Scandium 44.956
22 Ti Titanium 47.867
23 V Vanadium 50.942
24 Cr Chromium 51.996
25 Mn Manganese 54.938
26 Fe Iron 55.845
27 Co Cobalt 58.933
28 Ni Nickel 58.693
29 Cu Copper 63.546
30 Zn Zinc 65.38
31 Ga Gallium 69.723
32 Ge Germanium 72.630
33 As Arsenic 74.922
34 Se Selenium 78.971
35 Br Bromine 79.904
36 Kr Krypton 83.798
37 Rb Rubidium 85.468
38 Sr Strontium 87.62
39 Y Yttrium 88.906
40 Zr Zirconium 91.224
41 Nb Niobium 92.906
42 Mo Molybdenum 95.95
43 Tc Technetium 98
44 Ru Ruthenium 101.07
45 Rh Rhodium 102.91
46 Pd Palladium 106.42
47 Ag Silver 107.87
48 Cd Cadmium 112.41
49 In Indium 114.82
50 Sn Tin 118.71
51 Sb Antimony 121.76
52 Te Tellurium 127.60
53 I Iodine 126.90
54 Xe Xenon 131.29
55 Cs Caesium 132.91
56 Ba Barium 137.33
57 La Lanthanum 138.91
58 Ce Cerium 140.12
59 Pr Praseodymium 140.91
60 Nd Neodymium 144.24
61 Pm Promethium 145
62 Sm Samarium 150.36
63 Eu Europium 151.96
64 Gd Gadolinium 157.25
65 Tb Terbium 158.93
66 Dy Dysprosium 162.50
67 Ho Holmium 164.93
68 Er Erbium 167.26
69 Tm Thulium 168.93
70 Yb Ytterbium 173.05
71 Lu Lutetium 174.97
72 Hf Hafnium 178.49
73 Ta Tantalum 180.95
74 W Tungsten 183.84
75 Re Rhenium 186.21
76 Os Osmium 190.23
77 Ir Iridium 192.22
78 Pt Platinum 195.08
79 Au Gold 196.97
80 Hg Mercury 200.59
81 Tl Thallium 204.38
82 Pb Lead 207.2
83 Bi Bismuth 208.98
84 Po Polonium 209
85 At Astatine 210
86 Rn Radon 222
87 Fr Francium 223
88 Ra Radium 226
89 Ac Actinium 227
90 Th Thorium 232.04
91 Pa Protactinium 231.04
92 U Uranium 238.03
93 Np Neptunium 237
94 Pu Plutonium 244
95 Am Americium 243
96 Cm Curium 247
97 Bk Berkelium 247
98 Cf Californium 251
99 Es Einsteinium 252
100 Fm Fermium 257
101 Md Mendelevium 258
102 No Nobelium 259
103 Lr Lawrencium 266
104 Rf Rutherfordium 267
105 Db Dubnium 268
106 Sg Seaborgium 269
107 Bh Bohrium 270
108 Hs Hassium 269
109 Mt Meitnerium 278
110 Ds Darmstadtium 281
111 Rg Roentgenium 282
112 Cn Copernicium 285
113 Nh Nihonium 286
114 Fl Flerovium 289
115 Mc Moscovium 290
116 Lv Livermorium 293
117 Ts Tennessine 294
118 Og Oganesson 294`;

// [mass, abundance] pairs. Abundances are fractions summing to ~1.
const ISOTOPES = {
  H: [[1.00782503207, 0.999885], [2.0141017778, 0.000115]],
  Li: [[6.015122795, 0.0759], [7.01600455, 0.9241]],
  B: [[10.0129370, 0.199], [11.0093054, 0.801]],
  C: [[12.0, 0.9893], [13.0033548378, 0.0107]],
  N: [[14.0030740048, 0.99636], [15.0001088982, 0.00364]],
  O: [[15.99491461956, 0.99757], [16.99913170, 0.00038], [17.9991610, 0.00205]],
  F: [[18.99840322, 1]],
  Na: [[22.9897692809, 1]],
  Mg: [[23.985041700, 0.7899], [24.98583692, 0.1], [25.982592929, 0.1101]],
  Al: [[26.98153863, 1]],
  Si: [[27.9769265325, 0.92223], [28.976494700, 0.04685], [29.97377017, 0.03092]],
  P: [[30.97376163, 1]],
  S: [[31.97207100, 0.9499], [32.97145876, 0.0075], [33.96786690, 0.0425], [35.96708076, 0.0001]],
  Cl: [[34.96885268, 0.7576], [36.96590259, 0.2424]],
  K: [[38.96370668, 0.932581], [39.96399848, 0.000117], [40.96182576, 0.067302]],
  Ca: [[39.96259098, 0.96941], [41.95861801, 0.00647], [42.9587666, 0.00135], [43.9554818, 0.02086]],
  Fe: [[53.9396105, 0.05845], [55.9349375, 0.91754], [56.9353940, 0.02119], [57.9332756, 0.00282]],
  Cu: [[62.9295975, 0.6915], [64.9277895, 0.3085]],
  Zn: [[63.9291422, 0.4863], [65.9260334, 0.2790], [66.9271273, 0.0410], [67.9248442, 0.1875], [69.9253193, 0.0062]],
  Se: [[73.9224764, 0.0089], [75.9192136, 0.0937], [76.9199140, 0.0763], [77.9173091, 0.2377], [79.9165213, 0.4961], [81.9166994, 0.0873]],
  Br: [[78.9183371, 0.5069], [80.9162906, 0.4931]],
  Sn: [[115.901741, 0.1454], [116.902952, 0.0768], [117.901603, 0.2422], [118.903308, 0.0859], [119.9021947, 0.3258], [121.9034390, 0.0463], [123.9052739, 0.0579]],
  I: [[126.904473, 1]],
  Pt: [[193.9626803, 0.32967], [194.9647911, 0.33832], [195.9649515, 0.25242], [197.967893, 0.07163]],
  Pd: [[103.904036, 0.1114], [104.905085, 0.2233], [105.903486, 0.2733], [107.903892, 0.2646], [109.905153, 0.1172]],
};

// Readable atom label colors (light background). Others fall back to DEFAULT_COLOR.
const COLORS = {
  H: '#555555', C: '#222222', N: '#2244cc', O: '#dd2211', F: '#1e9e1e', Cl: '#1e9e1e',
  Br: '#a0331f', I: '#8a1a9a', S: '#b8960b', P: '#e07000', B: '#c76d6d', Si: '#8a6d4f',
  Li: '#9933cc', Na: '#7a3dc2', K: '#7a3dc2', Mg: '#3c8a2f', Fe: '#c45f1c', Se: '#b07600',
};
export const DEFAULT_COLOR = '#333333';

// Default valences (lowest first). Used for implicit hydrogen calculation.
export const VALENCES = {
  H: [1], B: [3], C: [4], N: [3, 5], O: [2], F: [1], Si: [4], P: [3, 5], S: [2, 4, 6],
  Cl: [1, 3, 5, 7], Ge: [4], As: [3, 5], Se: [2, 4, 6], Br: [1, 3, 5], Te: [2, 4, 6], I: [1, 3, 5, 7],
  Al: [3], Sn: [2, 4],
};

function tablePosition(z) {
  // Returns { row, col } (1-based) for a standard 18-column periodic table layout
  // with the f-block placed in rows 9 (lanthanides) and 10 (actinides).
  if (z === 1) return { row: 1, col: 1 };
  if (z === 2) return { row: 1, col: 18 };
  const blocks = [[3, 10, 2], [11, 18, 3]];
  for (const [s, e, row] of blocks) {
    if (z >= s && z <= e) {
      const i = z - s;
      return { row, col: i < 2 ? i + 1 : i + 11 };
    }
  }
  if (z <= 36) return { row: 4, col: z - 18 };
  if (z <= 54) return { row: 5, col: z - 36 };
  const heavy = (start, row) => {
    const i = z - start; // 0-based within period
    if (i < 2) return { row, col: i + 1 };
    if (i <= 16) return { row: row + 3, col: i + 1 }; // f-block row (La..Lu / Ac..Lr), cols 3..17
    return { row, col: i - 13 };
  };
  if (z <= 86) return heavy(55, 6);
  return heavy(87, 7);
}

export const ELEMENTS = {};
export const ELEMENT_LIST = [];
for (const line of RAW.split('\n')) {
  const [z, symbol, name, mass] = line.trim().split(/\s+/);
  const Z = Number(z);
  const iso = ISOTOPES[symbol] || null;
  let mono = Number(mass);
  if (iso) mono = iso.reduce((best, cur) => (cur[1] > best[1] ? cur : best))[0];
  const el = {
    z: Z, symbol, name, mass: Number(mass), monoisotopic: mono, isotopes: iso,
    color: COLORS[symbol] || DEFAULT_COLOR, valences: VALENCES[symbol] || null,
    ...tablePosition(Z),
  };
  ELEMENTS[symbol] = el;
  ELEMENT_LIST.push(el);
}

export function getElement(symbol) {
  return ELEMENTS[symbol] || null;
}

export function isElement(symbol) {
  return Object.prototype.hasOwnProperty.call(ELEMENTS, symbol);
}

export function elementColor(symbol) {
  return COLORS[symbol] || DEFAULT_COLOR;
}

const PAULING = {
  H: 2.20, Li: 0.98, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98, Na: 0.93, Mg: 1.31, Al: 1.61,
  Si: 1.90, P: 2.19, S: 2.58, Cl: 3.16, K: 0.82, Se: 2.55, Br: 2.96, I: 2.66, Sn: 1.96,
};
export function electronegativity(symbol) {
  return PAULING[symbol] ?? 1.8;
}

/** Number of valence electrons for main-group elements (used for aromaticity / radicals). */
export function valenceElectrons(symbol) {
  const el = ELEMENTS[symbol];
  if (!el) return 0;
  const c = el.col;
  if (c <= 2) return c;
  if (c >= 13) return c - 10;
  return 2;
}
