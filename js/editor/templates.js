// Ring templates, abbreviations (functional groups) and the built-in
// structure library. Library entries are SMILES so contributors can easily
// add more.

/** Ring tool definitions. `aromatic` rings get alternating double bonds. */
export const RING_TEMPLATES = [
  { id: 'ring3', size: 3, label: 'Cyclopropane' },
  { id: 'ring4', size: 4, label: 'Cyclobutane' },
  { id: 'ring5', size: 5, label: 'Cyclopentane' },
  { id: 'ring6', size: 6, label: 'Cyclohexane' },
  { id: 'ring7', size: 7, label: 'Cycloheptane' },
  { id: 'ring8', size: 8, label: 'Cyclooctane' },
  { id: 'benzene', size: 6, label: 'Benzene', aromatic: true },
  { id: 'cp', size: 5, label: 'Cyclopentadiene', aromatic: true },
];

/**
 * Abbreviations typed on an atom (Text/Label tool or hotkey). The first atom of
 * the SMILES is the attachment point. Typing e.g. "OMe" on a terminal atom
 * replaces it with the expanded group.
 */
export const ABBREVIATIONS = {
  Me: 'C', Et: 'CC', Pr: 'CCC', nPr: 'CCC', iPr: 'C(C)C', Bu: 'CCCC', nBu: 'CCCC', tBu: 'C(C)(C)C',
  iBu: 'CC(C)C', sBu: 'C(C)CC', Ph: 'c1ccccc1', Bn: 'Cc1ccccc1', Bz: 'C(=O)c1ccccc1', Ac: 'C(C)=O',
  OAc: 'OC(C)=O', OMe: 'OC', OEt: 'OCC', OPh: 'Oc1ccccc1', OBn: 'OCc1ccccc1', OH: 'O', SH: 'S', SMe: 'SC',
  NH2: 'N', NHMe: 'NC', NMe2: 'N(C)C', NEt2: 'N(CC)CC', NHAc: 'NC(C)=O', NO2: '[N+](=O)[O-]', CN: 'C#N',
  NC: '[N+]#[C-]', CHO: 'C=O', COOH: 'C(=O)O', CO2H: 'C(=O)O', COOMe: 'C(=O)OC', CO2Me: 'C(=O)OC',
  CO2Et: 'C(=O)OCC', COOEt: 'C(=O)OCC', CONH2: 'C(N)=O', COCl: 'C(Cl)=O', CF3: 'C(F)(F)F', CCl3: 'C(Cl)(Cl)Cl',
  SO3H: 'S(=O)(=O)O', SO2Me: 'S(C)(=O)=O', Ms: 'S(C)(=O)=O', OMs: 'OS(C)(=O)=O', Ts: 'S(=O)(=O)c1ccc(C)cc1',
  OTs: 'OS(=O)(=O)c1ccc(C)cc1', Tf: 'S(=O)(=O)C(F)(F)F', OTf: 'OS(=O)(=O)C(F)(F)F', Boc: 'C(=O)OC(C)(C)C',
  Cbz: 'C(=O)OCc1ccccc1', TMS: '[Si](C)(C)C', OTMS: 'O[Si](C)(C)C', TBS: '[Si](C)(C)C(C)(C)C',
  OTBS: 'O[Si](C)(C)C(C)(C)C', N3: 'N=[N+]=[N-]', Vinyl: 'C=C', Allyl: 'CC=C', CCH: 'C#C', Cy: 'C1CCCCC1',
  Bpin: 'B1OC(C)(C)C(C)(C)O1', PPh3: '[P+](c1ccccc1)(c1ccccc1)c1ccccc1',
};

/** Built-in library of common structures, grouped by category. */
export const STRUCTURE_LIBRARY = [
  {
    category: 'Common',
    items: [
      ['Benzene', 'c1ccccc1'], ['Toluene', 'Cc1ccccc1'], ['Phenol', 'Oc1ccccc1'], ['Aniline', 'Nc1ccccc1'],
      ['Benzoic acid', 'OC(=O)c1ccccc1'], ['Benzaldehyde', 'O=Cc1ccccc1'], ['Acetophenone', 'CC(=O)c1ccccc1'],
      ['Nitrobenzene', '[O-][N+](=O)c1ccccc1'], ['Styrene', 'C=Cc1ccccc1'], ['Ethanol', 'CCO'], ['Acetone', 'CC(C)=O'],
      ['Acetic acid', 'CC(=O)O'], ['Ethyl acetate', 'CCOC(C)=O'], ['Diethyl ether', 'CCOCC'], ['DMSO', 'CS(C)=O'],
      ['DMF', 'CN(C)C=O'], ['THF', 'C1CCOC1'], ['Dichloromethane', 'ClCCl'], ['Acetonitrile', 'CC#N'], ['Pyridine', 'c1ccncc1'],
    ],
  },
  {
    category: 'Rings & heterocycles',
    items: [
      ['Naphthalene', 'c1ccc2ccccc2c1'], ['Anthracene', 'c1ccc2cc3ccccc3cc2c1'], ['Phenanthrene', 'c1ccc2c(c1)ccc1ccccc12'],
      ['Indole', 'c1ccc2[nH]ccc2c1'], ['Quinoline', 'c1ccc2ncccc2c1'], ['Pyrrole', 'c1cc[nH]c1'], ['Furan', 'c1ccoc1'],
      ['Thiophene', 'c1ccsc1'], ['Imidazole', 'c1c[nH]cn1'], ['Pyrimidine', 'c1cncnc1'], ['Purine', 'c1ncc2[nH]cnc2n1'],
      ['Piperidine', 'C1CCNCC1'], ['Morpholine', 'C1COCCN1'], ['Piperazine', 'C1CNCCN1'], ['Cyclohexanone', 'O=C1CCCCC1'],
      ['Adamantane', 'C1C2CC3CC1CC(C2)C3'], ['Norbornane', 'C1CC2CCC1C2'], ['Steroid core', 'C1CCC2C(C1)CCC1C2CCC2CCCC21'],
    ],
  },
  {
    category: 'Drugs & natural products',
    items: [
      ['Aspirin', 'CC(=O)Oc1ccccc1C(=O)O'], ['Paracetamol', 'CC(=O)Nc1ccc(O)cc1'], ['Ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'],
      ['Caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'], ['Nicotine', 'CN1CCCC1c1cccnc1'], ['Vanillin', 'COc1cc(C=O)ccc1O'],
      ['Menthol', 'CC(C)C1CCC(C)CC1O'], ['Cholesterol', 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'],
      ['Penicillin G', 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O'], ['Morphine', 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O'],
      ['Dopamine', 'NCCc1ccc(O)c(O)c1'], ['Serotonin', 'NCCc1c[nH]c2ccc(O)cc12'], ['Adrenaline', 'CNCC(O)c1ccc(O)c(O)c1'],
      ['Glucose', 'OCC1OC(O)C(O)C(O)C1O'], ['Ascorbic acid', 'OCC(O)C1OC(=O)C(O)=C1O'],
    ],
  },
  {
    category: 'Amino acids',
    items: [
      ['Glycine', 'NCC(=O)O'], ['Alanine', 'CC(N)C(=O)O'], ['Valine', 'CC(C)C(N)C(=O)O'], ['Leucine', 'CC(C)CC(N)C(=O)O'],
      ['Isoleucine', 'CCC(C)C(N)C(=O)O'], ['Proline', 'OC(=O)C1CCCN1'], ['Phenylalanine', 'NC(Cc1ccccc1)C(=O)O'],
      ['Tyrosine', 'NC(Cc1ccc(O)cc1)C(=O)O'], ['Tryptophan', 'NC(Cc1c[nH]c2ccccc12)C(=O)O'], ['Serine', 'NC(CO)C(=O)O'],
      ['Threonine', 'CC(O)C(N)C(=O)O'], ['Cysteine', 'NC(CS)C(=O)O'], ['Methionine', 'CSCCC(N)C(=O)O'],
      ['Asparagine', 'NC(=O)CC(N)C(=O)O'], ['Glutamine', 'NC(=O)CCC(N)C(=O)O'], ['Aspartic acid', 'NC(CC(=O)O)C(=O)O'],
      ['Glutamic acid', 'NC(CCC(=O)O)C(=O)O'], ['Lysine', 'NCCCCC(N)C(=O)O'], ['Arginine', 'NC(=N)NCCCC(N)C(=O)O'],
      ['Histidine', 'NC(Cc1c[nH]cn1)C(=O)O'],
    ],
  },
  {
    category: 'Nucleobases & sugars',
    items: [
      ['Adenine', 'Nc1ncnc2[nH]cnc12'], ['Guanine', 'Nc1nc2[nH]cnc2c(=O)[nH]1'], ['Cytosine', 'Nc1cc[nH]c(=O)n1'],
      ['Thymine', 'Cc1c[nH]c(=O)[nH]c1=O'], ['Uracil', 'O=c1cc[nH]c(=O)[nH]1'], ['Ribose', 'OCC1OC(O)C(O)C1O'],
      ['Fructose', 'OCC1(O)OCC(O)C(O)C1O'],
    ],
  },
];

/** Offline name -> SMILES dictionary (used before falling back to PubChem). */
export function localNameLookup(name) {
  const key = name.trim().toLowerCase();
  for (const group of STRUCTURE_LIBRARY) {
    for (const [n, smi] of group.items) if (n.toLowerCase() === key) return smi;
  }
  const extra = {
    water: 'O', methane: 'C', ethane: 'CC', propane: 'CCC', butane: 'CCCC', pentane: 'CCCCC', hexane: 'CCCCCC',
    heptane: 'CCCCCCC', octane: 'CCCCCCCC', ethylene: 'C=C', ethene: 'C=C', acetylene: 'C#C', methanol: 'CO',
    propanol: 'CCCO', isopropanol: 'CC(C)O', formaldehyde: 'C=O', acetaldehyde: 'CC=O', 'formic acid': 'OC=O',
    ammonia: 'N', urea: 'NC(N)=O', chloroform: 'ClC(Cl)Cl', cyclohexane: 'C1CCCCC1', cyclopentane: 'C1CCCC1',
    'carbon dioxide': 'O=C=O', 'hydrogen peroxide': 'OO', glycerol: 'OCC(O)CO', 'lactic acid': 'CC(O)C(=O)O',
    'citric acid': 'OC(=O)CC(O)(CC(=O)O)C(=O)O', acetaminophen: 'CC(=O)Nc1ccc(O)cc1', xylene: 'Cc1ccccc1C',
    anisole: 'COc1ccccc1', benzonitrile: 'N#Cc1ccccc1', chlorobenzene: 'Clc1ccccc1', bromobenzene: 'Brc1ccccc1',
  };
  return extra[key] || null;
}
