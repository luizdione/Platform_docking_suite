# Docking Suite

An umbrella platform for structure-based antiparasitic screening: a declarative registry of
targets, ligands and protocols, four docking/simulation engines behind one interface, and a
dashboard that reads from it.

**This repository contains code and configuration only. It contains no research results.**
Every pose, trajectory, figure, log and report the platform produces is written to a data
tree *outside* the repository, pointed at by `DOCKING_DATA_DIR`. The two campaigns that
ship with it are examples built entirely from public PDB entries and published reference
compounds.

---

## What problem it solves

Most docking work starts as a script for one protein and one compound. It grows a second
target, then a human counterpart for selectivity, then a molecular dynamics step — and the
assumptions of the first system stay buried in the code. A force-field file written for a
flavoenzyme gets reused on a protein that has no flavin. A grid box silently falls back to
the centre of mass of an oligomer, which is solvent. Nothing errors out; the numbers just
stop meaning what you think they mean.

This platform's answer is that **the peculiarities of a system are configuration, not
code**. A target declares where its site is and which cofactors are structural. A ligand
declares whether it is a macrocycle, because a 19-membered ring treated as rigid gives a
wrong pose by construction. A **pair** — target × ligand, scoped to a campaign — declares
what only exists in the combination: the oligomeric state the site requires, whether a
measured affinity exists and in which state it was measured, whether the engine's score is
calibrated for this system at all.

Code reads the registry. Code contains no domain dictionary. `docs/REGISTRY.md` explains
why each field lives where it does.

## What it does

| Area | Module |
|---|---|
| Docking | AutoDock Vina (`engines/vina.py`), HADDOCK3 restraint-guided (`engines/haddock3.py`) |
| Molecular dynamics | OpenMM, with replicas and resume (`engines/openmm_md.py`, `md_simulation.py`, `replica_driver.py`) |
| Structure / affinity prediction | Boltz-2 (`engines/boltz2.py`) |
| ADMET | physicochemical descriptors, Lipinski/Veber, QED, permeability, CYP alerts (`admet_module.py`) |
| Toxicity and selectivity | in-silico alerts, GHS class, selectivity index from *your* assay data (`toxicity_module.py`) |
| Drug design | FBDD on a published R-group catalogue, BRICS hybridisation, genetic algorithm (`drug_design_module.py`) |
| Formulation | β-cyclodextrin inclusion feasibility (`cyclodex_module.py`) |
| Analysis | comparative figures and results table (`analysis.py`) |
| Provenance | plan/run/status/resume with a registry snapshot per run (`runner.py`) |
| Dashboard | Flask, port 5001 (`dashboard_v2.py`) |

---

## What ships with it

### Targets — 11 structures, all public PDB entries

No predicted model is included: every target declares `structure.source: pdb`, so nothing
here rests on a confidence score you would have to defend. The coordinate files are **not**
redistributed — the platform downloads them from the RCSB PDB on first use.

**Parasite enzymes**

| ID | PDB | Res. (Å) | Enzyme | Organism |
|---|---|---|---|---|
| `TcTR` | [1BZL](https://www.rcsb.org/structure/1BZL) | 2.40 | Trypanothione reductase | *Trypanosoma cruzi* |
| `TcCruzain` | [3KKU](https://www.rcsb.org/structure/3KKU) | 1.28 | Cruzain (cysteine protease) | *Trypanosoma cruzi* |
| `TcFeSOD` | [4DVH](https://www.rcsb.org/structure/4DVH) | 1.50 | Fe-superoxide dismutase | *Trypanosoma cruzi* |
| `LiTR` | [2JK6](https://www.rcsb.org/structure/2JK6) | 2.95 | Trypanothione reductase | *Leishmania infantum* |
| `LmDHODH` | [3C61](https://www.rcsb.org/structure/3C61) | 2.13 | Dihydroorotate dehydrogenase | *Leishmania donovani* |
| `LmPTR1` | [2QHX](https://www.rcsb.org/structure/2QHX) | 2.20 | Pteridine reductase 1 | *Leishmania major* |

**Human counterparts — the selectivity counter-targets**

| ID | PDB | Res. (Å) | Enzyme | Pairs with |
|---|---|---|---|---|
| `hGR` | [3DJJ](https://www.rcsb.org/structure/3DJJ) | 1.05 | Glutathione reductase | `TcTR`, `LiTR` |
| `hCatL` | [2XU3](https://www.rcsb.org/structure/2XU3) | 1.40 | Cathepsin L | `TcCruzain` |
| `hMnSOD` | [1LUV](https://www.rcsb.org/structure/1LUV) | 1.85 | Mn-superoxide dismutase | `TcFeSOD` |
| `hDHODH` | [1D3H](https://www.rcsb.org/structure/1D3H) | 1.60 | Dihydroorotate dehydrogenase | `LmDHODH` |

**Specificity counter-target**

| ID | PDB | Res. (Å) | Assembly |
|---|---|---|---|
| `tubulin_ab` | [1SA0](https://www.rcsb.org/structure/1SA0) | 3.58 | α/β tubulin heterodimer, colchicine site |

The parasite/human pairing is the point, not decoration. A docking score alone says nothing
about selectivity; the platform docks the parasite enzyme and its closest human counterpart
under the same protocol and reads selectivity from the difference.

### Reference ligands — 7 published compounds

| ID | Compound | Class | Reference structure |
|---|---|---|---|
| `beta-lapachone` | β-Lapachone | naphthoquinone, natural product | — |
| `plumbagin` | Plumbagin | naphthoquinone, natural product | — |
| `colchicina` | Colchicine | tubulin, colchicine site | [1SA0](https://www.rcsb.org/structure/1SA0) (DAMA-colchicine, `CN2`) |
| `citocalasina-d` | Cytochalasin D | actin-binding macrolide | [3EKS](https://www.rcsb.org/structure/3EKS) (`CY9`) |
| `jasplaquinolida` | Jasplakinolide | actin-binding macrocycle | [7PM3](https://www.rcsb.org/structure/7PM3) (`9UE`) |
| `latrunculina-a` | Latrunculin A | actin-binding macrolide | [1ESV](https://www.rcsb.org/structure/1ESV) (`LAR`) |
| `latrunculina-b` | Latrunculin B | actin-binding macrolide | [2Q0U](https://www.rcsb.org/structure/2Q0U) (`LAB`) |

`populate_ligand_db.py` seeds a further 18 published compounds — clinical reference drugs
(benznidazole, nifurtimox, miltefosine, amphotericin B, atovaquone), natural quinones and
described synthetic derivatives — each with its primary citation. Comparators with
provenance, so a score has something to be compared against.

### Example campaigns

- **`DEMO01`** — antiparasitic selectivity screen: six parasite enzymes against four human
  counterparts, with β-lapachone and plumbagin as ligands.
- **`DEMO02`** — specificity cross-docking at the tubulin colchicine site. Colchicine is
  the **positive** control and must reproduce its own crystallographic pose; the four
  actin-binding natural products are **negative** controls, because none has a site there.
  Both arms are needed: without the positive, a weak score on the negatives might just mean
  the box is wrong.

---

## Install

RDKit and Open Babel are the reason conda is the recommended route; pip alone is fragile
for them on Windows.

```bash
conda create -n docking python=3.11
conda activate docking
conda install -c conda-forge numpy scipy matplotlib pandas scikit-learn biopython rdkit openbabel
pip install -r requirements.txt
```

On Windows, `setup_bioinfo_env.bat` does the second half for an already-active environment.
It reads `CONDA_ROOT`/`ENV_NAME` from the environment and contains no hard-coded paths.

### External engines

None of these is bundled; each keeps its own licence and citation requirement.

| Engine | Needed for | Notes |
|---|---|---|
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) ≥ 1.2 | docking | the one hard requirement for a docking run |
| [OpenMM](https://openmm.org) ≥ 8.4 | molecular dynamics | NVIDIA GPU strongly recommended; CPU works and is slow |
| [HADDOCK3](https://github.com/haddocking/haddock3) | restraint-guided docking | academic registration required |
| [Boltz-2](https://github.com/jwohlwend/boltz) | structure/affinity prediction | optional |

Put the binaries on `PATH`, or register the paths in `registry/toolchain.json` (it ships
empty, with an example). `python engines/toolchain.py` reports what it found.

### Where results go — set this first

```bash
export DOCKING_DATA_DIR=/path/to/your/data        # Linux/macOS
set DOCKING_DATA_DIR=D:\data\docking              # Windows
```

With the variable unset, the root is `./docking_data` relative to the working directory,
and `.gitignore` blocks it. The layout, per campaign:

```
$DOCKING_DATA_DIR/<CAMPAIGN>/
    proteins/                       prepared receptors
    ligand/<ligand>/                prepared ligands
    runs/<run_id>/run.json          provenance: registry snapshot + hashes
    results/<ligand>/results.json   aggregated result
    figures/  logs/  reports/
```

`python paths.py` prints where the current campaign reads and writes. It is the first thing
to run when the dashboard says there is no result.

Other optional variables: `DOCKING_CAMPAIGN` (override the current campaign),
`DOCKING_PYTHON` (interpreter for sub-runs), `OBABEL_EXE`, `NVIDIA_SMI`.

---

## Run the example campaign

Check the registry and see the resolved plan without executing anything — `plan` is pure,
it reads and decides and touches no file:

```bash
python registry.py                                    # validate; [ok] registry valido
python runner.py plan --campaign DEMO01 --kind docking
```

Then run it. Vina must be reachable:

```bash
python runner.py run --campaign DEMO01 --kind docking
```

Each run writes `runs/<run_id>/run.json` with the full registry snapshot and its SHA-256,
the resolved toolchain, and per task the command line, input/output hashes, wall time and
every quality gate with its verdict. A failed gate gives `status: gate_failed`, not an
exception — the run records that it was refused and why. `runner.py resume` picks up from
that manifest.

The dashboard, for inspection and figures:

```bash
python dashboard_v2.py          # http://127.0.0.1:5001
```

Tests:

```bash
python -m pytest tests/ -q      # 63 pass; 3 skip, each stating what it needs
```

The three skips are integration gates: they need prepared receptors, a geometric golden
(`python tools/make_grid_golden.py`) or a second project's data. They say so rather than
passing without having measured anything.

---

## Container image

```bash
docker pull ghcr.io/luizdione/platform_docking_suite:latest

# one job, with your data mounted
docker run --gpus all -v /path/to/your/data:/data \
  ghcr.io/luizdione/platform_docking_suite:latest \
  runner.py run --campaign DEMO01 --kind docking

# dashboard
docker run -p 5001:5001 -v /path/to/your/data:/data \
  ghcr.io/luizdione/platform_docking_suite:latest
```

The image carries code and dependencies, plus AutoDock Vina. It carries **no data**:
`.dockerignore` blocks every data path, and `/data` is an empty mount point, so anything
written there without a volume dies with the container. `DOCKING_DATA_DIR` defaults to
`/data` inside the image.

### How the image is published

By GitHub Actions — `.github/workflows/publicar-imagem.yml`, which runs on the GitHub
runners rather than on a workstation.

The tagging rule exists so that a version tag never lies about its contents:

| Event | Tags applied |
|---|---|
| push to `main` touching the Dockerfile or dependencies | `main`, `sha-<short>` |
| **published release** | the version from the tag, and `latest` |
| manual run | **only** `sha-<short>` of the ref you ask for |

A manual run deliberately applies no version and does not move `latest`: `checkout`
defaults to the default branch, so a manual run asking for a version tag would rebuild
that tag from `main` and the tag would then serve content that is not that version.

Before anything is pushed, the job **checks the image**: it builds locally, exports the
filesystem, and fails if any structure or trajectory file (`.pdb`, `.pdbqt`, `.dcd`,
`.cif`, `.npz`, …) or any data directory appears under `/app`; then it runs
`vina --version` and `registry.py` inside the container. Only then does it push.
`.dockerignore` should already have filtered those paths — the gate exists because an
image in a registry has already been distributed.

To publish a version, create a release for the tag; the workflow does the rest.

---

## Documentation

| Document | What is in it |
|---|---|
| [`docs/PROTOCOL.md`](docs/PROTOCOL.md) | the scientific protocol: target table, docking and MD parameters, quality gates, how to read the metrics, and the traps that produce plausible wrong numbers |
| [`docs/REGISTRY.md`](docs/REGISTRY.md) | the declarative registry: hierarchy, why the pair is where the peculiarity lives, which fields change pipeline behaviour |
| [`docs/DRUG_DESIGN_MANUAL.md`](docs/DRUG_DESIGN_MANUAL.md) | the drug-design module: published R-group catalogue, BRICS, genetic algorithm, the empirical QSAR model and its limits |

---

## Licence

Custom licence — see [`LICENSE`](LICENSE) for the terms, which prevail over this summary.

- **Personal use and teaching: free**, no permission needed.
- **Research use requires the author's prior written authorisation**, by e-mail. That
  includes any use that generates or analyses data for a dissertation, thesis, report,
  grant application, preprint, article, poster or patent. Running the platform only to
  evaluate it does not.
- **Commercial use: not permitted** under this licence. Contact the author to negotiate.
- Redistribution, modification and derivative works: not permitted without written
  authorisation.

Third-party components keep their own licences — see
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). 3Dmol.js is bundled (BSD-3); the
docking, MD and cheminformatics engines are called, not distributed.

Your data stays yours: the platform runs in your environment, and the author neither
collects nor accesses what it produces.

## Citation

If you use the platform in authorized research, cite it as in
[`CITATION.cff`](CITATION.cff), and cite the engines you actually ran — Vina, OpenMM,
HADDOCK3, Boltz-2, RDKit and Open Babel each have their own reference.

## Author

**Luiz Dione Barbosa De-Melo** (De-Melo, L.D.B.)
ORCID [0000-0003-2445-9943](https://orcid.org/0000-0003-2445-9943) ·
<luiz.melo@ifrj.edu.br>
Instituto Federal de Educação, Ciência e Tecnologia do Rio de Janeiro (IFRJ), Brazil
