# Third-party notices

The Docking Suite is distributed under its own licence (see `LICENSE`). That licence does
**not** apply to the third-party components listed here. Each one keeps its own licence,
which prevails over `LICENSE` for that component.

Two kinds of component are listed: code **bundled** in this repository, and external
programs the platform **calls** but does not ship.

---

## Bundled in this repository

### 3Dmol.js — `static/3Dmol-min.js`

A WebGL molecular viewer, bundled so the dashboard renders structures without a network
round-trip to a CDN.

- Upstream: <https://3dmol.org> · <https://github.com/3dmol/3Dmol.js>
- Licence: BSD 3-Clause
- Reference: Rego N., Koes D. *3Dmol.js: molecular visualization with WebGL.*
  Bioinformatics 2015;31(8):1322–1324. <https://doi.org/10.1093/bioinformatics/btu829>
- The minified bundle carries the pointer `/*! For license information please see
  3Dmol-min.js.LICENSE.txt */`; the full licence text of 3Dmol.js and of its own
  dependencies is published with the upstream build.

This file is redistributed unmodified. Do not edit it: to change versions, replace it with
a fresh upstream build.

---

## Loaded at runtime, not bundled

### Google Fonts — Inter, JetBrains Mono

The dashboard references the Google Fonts CSS API (`fonts.googleapis.com`). The font files
are not in this repository; they are fetched by the viewer's browser. Inter and
JetBrains Mono are both published under the SIL Open Font License 1.1.

Consequence worth knowing: with the dashboard open, the browser contacts
`fonts.googleapis.com`. On an isolated network the fonts simply fall back to the system
stack and the interface keeps working.

### IFRJ logo

The dashboard header links the institutional logo from `www.ifrj.edu.br`. It is an
institutional mark of the Instituto Federal de Educação, Ciência e Tecnologia do Rio de
Janeiro, used here to identify the institution of origin. It is not covered by `LICENSE`
and is not granted for reuse.

---

## External programs the platform calls

None of these is distributed here. Install them yourself; each keeps its own licence and
its own citation requirement.

| Program | Role in the platform | Licence | Where |
|---|---|---|---|
| AutoDock Vina | docking engine (`engines/vina.py`) | Apache 2.0 | <https://github.com/ccsb-scripps/AutoDock-Vina> |
| OpenMM | molecular dynamics (`engines/openmm_md.py`, `md_simulation.py`) | MIT / LGPL | <https://openmm.org> |
| HADDOCK3 | restraint-guided docking (`engines/haddock3.py`) | Apache 2.0, **academic registration required** | <https://github.com/haddocking/haddock3> |
| Boltz-2 | structure/affinity prediction (`engines/boltz2.py`) | MIT | <https://github.com/jwohlwend/boltz> |
| RDKit | cheminformatics (ADMET, toxicity, drug design) | BSD 3-Clause | <https://www.rdkit.org> |
| Open Babel | 3D generation and PDBQT preparation | GPL 2.0 | <https://openbabel.org> |
| Meeko | PDBQT preparation for Vina | LGPL 2.1 | <https://github.com/forlilab/Meeko> |
| Biopython | PDB parsing and structure handling | Biopython licence (MIT-like) | <https://biopython.org> |
| NumPy, SciPy, Matplotlib, pandas | numerics and figures | BSD 3-Clause | — |
| Flask | dashboard web server | BSD 3-Clause | — |

Open Babel is GPL 2.0. The platform invokes it as a **separate executable** through the
command line (`compound_db._find_obabel`), which is why the Docking Suite's own licence is
not affected by it. If you link against it, that is your obligation to assess, not a
permission granted here.

---

## Structural and chemical data

The target and ligand entries in `registry/` describe **public** data: PDB entries and
compounds described in the published literature, cited in the `provenance` and `notes`
fields of each file. The structure files themselves are not redistributed here — the
platform downloads them from the RCSB PDB when a campaign runs.

Terms of use of the source databases apply to that data and are theirs, not ours:

- RCSB Protein Data Bank — <https://www.rcsb.org/pages/usage-policy> (CC0 1.0)
- ChEMBL — <https://chembl.gitbook.io/chembl-interface-documentation/about> (CC BY-SA 3.0)
