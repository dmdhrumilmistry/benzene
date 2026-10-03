# Contributing to Benzene

Thanks for your interest in improving Benzene!

## Principles

- **Static and dependency-free.** The app must keep working as plain files on GitHub Pages: no bundler, no npm runtime dependencies, no backend.
- **Privacy first.** Never send user structures anywhere unless the user explicitly triggers it, as with the PubChem lookups today.
- **Chemistry stays DOM-free.** Everything in `js/core/` must run in Node so it can be unit tested.

## Development

```bash
python3 -m http.server 8080   # serve the app at http://localhost:8080
npm test                      # run unit tests (Node ≥ 20)
```

## Code style

- Modern JavaScript (ES2022), ES modules, 2-space indent, single quotes, semicolons.
- Add short JSDoc comments to exported functions.
- Keep data tables (NMR increments, templates, abbreviations) readable and commented. Many contributors will only ever touch those.

## Conventions in the molecule model

- Coordinates are in **bond-length units** (a standard bond is 1.0), with the **y axis pointing down**.
- Bonds are stored in Kekulé form (orders 1, 2, 3). Aromaticity is derived on demand by `js/core/rings.js`.
- `atom.hCount === null` means implicit hydrogens are computed from the valence.

## Pull requests

1. Fork the repository and create a feature branch.
2. Add or update tests in `tests/` for any change to `js/core/`.
3. Make sure `npm test` passes.
4. For UI changes, include a screenshot in the PR description.

## Ideas wanted

- E/Z and tetrahedral stereo perception (CIP labels), and stereo SMILES output
- CDXML import/export
- InChI generation (likely via an optional WebAssembly module)
- More accurate NMR increments, 2D NMR (COSY/HSQC) sketches, and IR band hints
- Reaction schemes: automatic layout of reactants → products
- Touch gestures (pinch zoom) and accessibility improvements
