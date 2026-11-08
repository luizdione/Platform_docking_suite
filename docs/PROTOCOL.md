# Scientific protocol

What the platform does to a protein–ligand pair, and why each step is the way it is. This
is the document to read before trusting a number the dashboard prints.

Two things it is **not**: a record of results (there are none in this repository), and a
source of truth about targets. Every target, ligand, pair and protocol lives in
`registry/` — see [REGISTRY.md](REGISTRY.md). The tables below were **generated from the
registry**, not transcribed into prose. That matters: an earlier version of this document
was maintained by hand and drifted, at one point attributing to a parasite cysteine
protease a PDB code that belongs to human glutathione reductase. Regenerate the tables
instead of editing them:

```bash
python registry.py --campaign DEMO01
```

---

## 1. Targets included

Eleven structures ship with the platform, all of them public PDB entries. No predicted
model is included: every target declares `structure.source: pdb`, so nothing here depends
on a confidence score you would have to argue about.

**Parasite enzymes**

| ID | PDB | Res. (Å) | Enzyme | Organism | Cofactors kept | Site definition | Box (Å) |
|---|---|---|---|---|---|---|---|
| `TcTR` | `1BZL` | 2.40 | Trypanothione reductase | *Trypanosoma cruzi* | FAD | `ligand_center` | 25×25×25 |
| `TcCruzain` | `3KKU` | 1.28 | Cruzain (cysteine protease) | *Trypanosoma cruzi* | — | `ligand_center` | 22×22×22 |
| `TcFeSOD` | `4DVH` | 1.50 | Fe-superoxide dismutase | *Trypanosoma cruzi* | FE, FE2 | `metal_center` | 22×22×22 |
| `LiTR` | `2JK6` | 2.95 | Trypanothione reductase | *Leishmania infantum* | FAD | `ligand_center` | 25×25×25 |
| `LmDHODH` | `3C61` | 2.13 | Dihydroorotate dehydrogenase | *Leishmania donovani* | FMN, FAD | `ligand_center` | 25×25×25 |
| `LmPTR1` | `2QHX` | 2.20 | Pteridine reductase 1 | *Leishmania major* | NAP, NDP, FAD | `ligand_center` | 22×22×22 |

**Human counterparts — the selectivity counter-targets**

| ID | PDB | Res. (Å) | Enzyme | Cofactors kept | Site definition | Box (Å) |
|---|---|---|---|---|---|---|
| `hGR` | `3DJJ` | 1.05 | Glutathione reductase | FAD | `ligand_center` | 25×25×25 |
| `hCatL` | `2XU3` | 1.40 | Cathepsin L | — | `residue_center` | 22×22×22 |
| `hMnSOD` | `1LUV` | 1.85 | Mn-superoxide dismutase | MN | `metal_center` | 22×22×22 |
| `hDHODH` | `1D3H` | 1.60 | Dihydroorotate dehydrogenase | FMN | `ligand_center` | 25×25×25 |

**Specificity counter-target**

| ID | PDB | Res. (Å) | Assembly | Site definition | Box (Å) |
|---|---|---|---|---|---|
| `tubulin_ab` | `1SA0` | 3.58 | α/β heterodimer | `explicit_center` | 27×25.6×22.5 |

### Why pair each parasite enzyme with a human one

A docking score on its own says nothing about selectivity. The platform's central pattern
is to dock the parasite enzyme **and its closest human counterpart under the same
protocol**, and to read selectivity from the difference. The example campaign declares
these pairs:

| Parasite target | Human counterpart | Why this is the right comparison |
|---|---|---|
| `TcTR` | `hGR` | Functionally equivalent flavoenzymes on different thiols — trypanothione vs glutathione. Trypanothione reductase has the wider active site, because it has to accommodate the spermidine bridge of its substrate. The pathway is absent in mammals. |
| `LiTR` | `hGR` | Same enzyme family in *Leishmania*; lets one run compare *T. cruzi* and *Leishmania* side by side. |
| `TcFeSOD` | `hMnSOD` | Iron- vs manganese-dependent superoxide dismutase. The Fe form is absent in humans, who use Cu/Zn and Mn enzymes. |
| `TcCruzain` | `hCatL` | Homologous cysteine proteases — the classic cross-reactivity risk for cruzain inhibitors. |
| `LmDHODH` | `hDHODH` | Same reaction, structurally distinct enzymes (class 1A vs class 2). |

> **On the `LiTR` proxy.** There is no crystal structure of trypanothione reductase from
> *L. amazonensis*. `2JK6` (*L. infantum*) is the standard proxy in the literature at ~95%
> sequence identity. That substitution is a modelling decision, not a fact, and it is
> recorded in the target's own `provenance` field so a reviewer can find it without
> reading the code. References: Baiocco et al., *Structure* 2009;17:1547; Camargo et al.,
> *Sci. Rep.* 2023;13:6887.

---

## 2. Docking

### 2.1 Ligand preparation

| Step | Setting |
|---|---|
| 3D generation | Open Babel (pybel), SMILES → 3D |
| Force field | MMFF94 |
| Optimisation | 1000 steps, conjugate gradients |
| Hydrogens | added |
| Output | PDBQT (via `compound_db.py`; falls back to RDKit + Meeko) |

The PDBQT is cached per ligand under `<DOCKING_DATA_DIR>/<CAMPAIGN>/ligand/<name>/`. A
ligand is prepared once per campaign, not once per target.

**A trap worth naming.** Open Babel writes PDBQT in *ligand* flavour, with `ROOT`/`BRANCH`
records. Fed to Vina as a receptor, that file is accepted and produces nonsense.
`docking_v2.py` strips those records when preparing a receptor, and `tests/` has a gate for
corrupt receptor PDBQTs left over by older pipelines. If you prepare a receptor by hand,
use UCSF Chimera → DockPrep → Save PDBQT.

### 2.2 Receptor preparation — what is removed and what must stay

**Removed:** water (HOH, WAT, DOD), buffer and cryoprotectant ions and molecules (SO4,
PO4, CL, GOL, EDO, MPD, …), and the co-crystallised ligand — which is read first, to place
the box, and only then deleted.

**Kept, and this is the part that silently ruins runs when it is wrong:**

| Kept | Why removing it breaks the result |
|---|---|
| FAD | The catalytic site of a trypanothione/glutathione reductase forms *on the flavin face*. Without FAD there is no pocket, only a cavity. |
| FMN | Essential cofactor of dihydroorotate dehydrogenase. |
| NADPH (NAP/NDP) | Part of the binding surface, not an accessory. |
| Fe²⁺, Mn²⁺, Zn²⁺, Mg²⁺, Ca²⁺ | In a metalloenzyme the metal *is* the active site. |

The default lists live in `registry/program.json` under `chemistry_defaults`, and any
target may override them in `chemistry.essential_cofactors`. The defaults are deliberately
biased toward flavoenzymes, because that is where the platform started — keeping that as a
*default* is correct; assuming it is universal is not. A target whose cofactor is
transferred from another structure rather than crystallised with it must declare
`cofactor_source: transferred` plus the donor and the alignment scope, and the registry
refuses the file otherwise.

### 2.3 Where the box goes

| Priority | Strategy | When it applies |
|---|---|---|
| 1 | `ligand_center` — centroid of the co-crystallised ligand | whenever the PDB entry has one |
| 2 | `metal_center` — the catalytic metal | Fe-SOD, Mn-SOD |
| 3 | `residue_center` — centroid of named active-site residues | apo structures |
| 4 | `derived_from_reference_complex` | no co-crystal of its own; the site is inferred by superposition onto a reference complex |
| 5 | `explicit_center` | coordinates recorded directly, for sites that none of the above reaches — e.g. an interface site |

Box sizes follow the site, not the target: ~25 Å per edge for the wide reductase site that
accommodates a bulky substrate, ~22 Å for a compact protease or metal site.

**Never let the box fall back to the centre of mass.** For an oligomeric receptor the
centre of mass is usually solvent, the docking still "works", and the result is garbage
with no error. `tests/test_oligomer_site_safety.py` enforces that every multi-protomer
target has an explicit `box_center_A`.

### 2.4 AutoDock Vina

```
exhaustiveness : 32   (production)  ·  8 (quick check)  ·  64 (maximum)
num_modes      : 10
energy_range   : 3 kcal/mol
scoring        : vina
```

Declared in `registry/protocols/vina_exh32.json`, and a pair may override it.

### 2.5 Restraint-guided docking (HADDOCK3)

Vina resolves an enclosed cavity well. It is the wrong tool when the site is defined by
*which residues face each other* rather than by a pocket: a protomer–protomer junction in a
filament, a dimer interface, or a site known from catalytic residues with no co-crystal.
For those, `registry/protocols/haddock3_default.json` drives HADDOCK3 with ambiguous
interaction restraints anchored on the declared residues.

---

## 3. Molecular dynamics

MD is not a better docking. It answers one question docking cannot: **does the pose
survive?** A pose that drifts out of the pocket in 5 ns was a scoring artefact.

### 3.1 Parameters

| Parameter | Value (`md_simulation.py`) |
|---|---|
| Force field | AMBER ff14SB + TIP3P |
| Solvent | TIP3P, 1 nm padding |
| Ionic strength | 150 mM NaCl, neutralised |
| Temperature | 300 K |
| Pressure | 1 bar, NPT, Monte Carlo barostat |
| Timestep | 4 fs with hydrogen-mass repartitioning (2 fs without) |
| Equilibration | NVT 50 ps → NPT 500 ps, restraints released gradually |
| Trajectory frequency | every 5 000 steps (~10 ps/frame) |
| Checkpoint | every 50 000 steps (~100 ps) |
| Minimum for a verdict | 5 ns |
| Production | 20–100 ns |

Two protocols ship, and they are **not** interchangeable:

- `md_ff14sb_hmr_300k` — ff14SB, HMR, 4 fs, 300 K. The default.
- `md_charmm36m_2fs_310k` — CHARMM36m, 2 fs, 310 K. Required when the target needs a
  **nucleotide cofactor** (ATP/ADP/Mg): ff14SB has no native templates for those residues,
  and a run that silently drops the nucleotide has lost the site.

### 3.2 Replicas: why n=1 is not a result

A single trajectory reports the luck of one velocity draw. `replica_driver.py` runs
independent replicas of the same starting pose — new tag per replica, so velocities are
reassigned in the warm-up and the Langevin integrator diverges — with per-run resume and an
idempotent state file:

```bash
python replica_driver.py --target TcTR --replicas 3
python replica_driver.py --jobs minhas_replicas.json    # full control
```

Each replica keeps a **distinct tag**. The script writes `md_<tag>.mode` and refuses to
resume a checkpoint from a different configuration, with an explicit error — mixing
checkpoints across configurations is the one failure that produces a plausible trajectory
of nothing.

### 3.3 Reading the trajectory

Use a library that reimages periodic boundaries. Measuring a distance in a viewer that does
not gives you the wrapping artefact, not the distance — typically a nonsense value in the
hundreds of ångström for two atoms that are in contact.

```python
import mdtraj as md
t = md.load("md_<tag>.dcd", top="minimized_<tag>.pdb")
d = md.compute_distances(t, pairs, periodic=True)   # minimum image
```

| Metric | Reading |
|---|---|
| Backbone Cα RMSD | ≤ 2.0 Å — the protein is stable; the run is interpretable |
| Ligand RMSD in the pocket | ≤ 2.0 Å — the pose held |
| Ligand → cofactor distance | mechanistic, per hypothesis; always with `periodic=True` |

`analyze_dcd_numpy.py` computes RMSD/RMSF without mdtraj when it is unavailable, including
the case where the topology has fewer atoms than the trajectory (different water counts
between builds): it falls back to the common prefix, which is deterministic for
protein + cofactor + ligand + ions.

### 3.4 GPU and power draw

MD saturates the GPU for hours. On a workstation fed by a modest UPS, a CPU-heavy job at
100% plus MD at 100% sums the peaks and can drop the host mid-production. Two mitigations,
in order: **serialise** the heavy jobs, and cap the GPU power limit
(`nvidia-smi -pl <W>`, requires administrator — the call fails by exit code, not by
exception, so a `try/except` around it reports success while doing nothing).
`run_md_safe.ps1` does both, plus the force-field gate below, and blocks sleep/hibernate
for the duration.

### 3.5 The force-field gate

`run_md_safe.ps1` runs `validate_ff.py` on the cofactor and ligand templates **before**
starting, and aborts if either fails. A missing or malformed template does not stop OpenMM:
it produces a system with the wrong parameters and runs for fourteen hours.

---

## 4. ADMET

`admet_module.py` computes physicochemical descriptors and classical empirical predictions
from the ligand SMILES with RDKit: Lipinski and Veber rules, QED, TPSA, estimated Caco-2
permeability, GI absorption, blood–brain barrier penetration, CYP structural alerts.

Two cautions that belong in the methods section of anything you write from this:

- These are **predictions from empirical models**, not measurements. They are a filter for
  prioritisation, not evidence.
- For an antiparasitic, "does not cross the blood–brain barrier" is usually the desirable
  answer, not a defect — the opposite of the reading in a CNS programme.

---

## 5. Toxicity and selectivity

`toxicity_module.py` computes in-silico alerts and, **when you supply your own
experimental values**, recomputes the selectivity index from them. The form ships empty:
the platform distributes no one's assay data.

### 5.1 Selectivity index

```
SI = CC50 (host cell) / IC50 (parasite)
```

| SI | Reading |
|---|---|
| > 50 | high selectivity |
| > 10 | favourable |
| 3–10 | moderate — a limitation to discuss, not to hide |
| < 3 | non-selective cytotoxicity |

An SI computed from a CC50 and an IC50 **measured in different assay formats** is not
comparable. Keep the pair consistent or say which is which.

### 5.2 Which CC50

| Line | Why | Priority |
|---|---|---|
| VERO (monkey kidney) | the reference in the antiparasitic literature | high |
| THP-1 (human macrophage) | the parasite lives *inside* the macrophage | high |
| J774.A1 (murine macrophage) | murine equivalent of THP-1 | medium |
| HepG2 (hepatocyte) | hepatic toxicity for an oral candidate | medium |

### 5.3 Which parasite IC50

| Form | Clinical relevance | Note |
|---|---|---|
| Intracellular trypomastigotes / amastigotes | highest | the form that actually matters |
| Axenic epimastigotes / promastigotes | screening | vector-stage form; **overestimates** activity |

Reporting an axenic IC50 as if it were the clinical one is the most common way a selectivity
index ends up flattering a compound.

### 5.4 GHS classification

| Class | Oral LD50 (rat) | Label |
|---|---|---|
| 1–2 | ≤ 50 mg/kg | fatal |
| 3 | 50–300 mg/kg | toxic |
| 4 | 300–2 000 mg/kg | harmful |
| 5–6 | > 2 000 mg/kg | little or no hazard |

### 5.5 Cyclodextrin inclusion

Poorly soluble quinones are classic candidates for inclusion complexation.
`cyclodex_module.py` estimates geometric fit against the β-cyclodextrin cavity (inner
diameter ~6.2 Å, depth ~7.8 Å) and the thermodynamics of inclusion, and reports the
stoichiometry it assumed.

Without RDKit the module **fails loudly** rather than falling back to assumed dimensions: a
fit verdict computed from a guessed molecular size looks exactly like a real one. β-CD is
FDA GRAS for oral use; any predicted stability constant still has to be confirmed by phase
solubility, DSC/FTIR and ¹H NMR before it means anything.

---

## 6. Provenance

Every run through `runner.py` writes `runs/<run_id>/run.json` into the campaign's data
tree, containing:

- the **full registry snapshot** used, plus its SHA-256 — so you can prove which
  configuration produced a result, and detect that it changed;
- the resolved toolchain (which Vina, which Python, which interpreter path);
- per task: command line, input and output hashes, wall time, and each quality gate with
  its verdict.

A failed gate yields `status: gate_failed`, not an exception — the run records that it was
refused and why. `runner.py resume` skips tasks already completed, keyed on that manifest.

Run `python runner.py plan --campaign DEMO01 --kind docking` to see the resolved plan
without executing anything. The plan is pure: it reads the registry and decides, and it
touches no file.

---

## 7. Quality gates

These are platform rules, not campaign preferences. They are in `runner.py` and a pair may
override individual thresholds through the registry.

| Gate | Threshold | What it prevents |
|---|---|---|
| `torsions_min` | ≥ 15 rotatable bonds in the ligand PDBQT | a PDBQT written in the wrong flavour, or a preparation that silently rigidified the molecule. **Rigid macrocycles are legitimately below it** — a 19-membered ring treated as rigid gives a wrong pose by construction, which is why `topology.is_macrocycle` forces `prep.pdbqt.rigid_macrocycles: false`, and those pairs set the threshold themselves. |
| `redock_rmsd` | ≤ 2.0 Å | the box or the preparation is wrong. If redocking the crystallographic ligand does not reproduce its own pose, nothing else in that target is worth reading. |
| `site_plddt` | ≥ 70 | a predicted model whose *site* is low confidence, even when the global score looks fine. |
| `protonation` | post-check required | protonation fails silently and the docking returns a wrong pose with no warning. |

RMSD is computed with `CalcRMS`, never `GetBestRMS`: the latter aligns before measuring and
therefore hides a wrong pose. A real 7.86 Å error once read as 0.40 Å that way, and the
registry now rejects any pair that asks for it.

---

## 8. Controls

A specificity matrix needs both arms. The example campaign `DEMO02` is built to show this:
colchicine docked into the colchicine site of tubulin is the **positive** arm and must
reproduce the pose from co-crystal `1SA0`; actin-binding natural products — cytochalasin D,
jasplakinolide, latrunculin A and B — docked into that same box are the **negative** arms,
because none of them has a site there.

Without the positive arm, a weak score on the negatives proves nothing: the box might
simply be wrong. Without the negatives, a strong score on the positive might be a scoring
function inflating with ligand volume. The registry enforces that any pair marked as a
control declares `type`, a `expected_binding` boolean consistent with that type, and a
`rationale` — and `tests/test_control_schema.py` is the gate.

---

## 9. Where results go

Nowhere inside this repository. Everything is written under `$DOCKING_DATA_DIR`:

```
<DOCKING_DATA_DIR>/<CAMPAIGN>/
    proteins/                       prepared receptors
    ligand/<ligand>/                prepared ligands
    runs/<run_id>/run.json          provenance
    results/<ligand>/results.json   aggregated result
    figures/<ligand>/<session>/
    logs/<ligand>/<session>/
    reports/<ligand>/<session>/     generated text
```

With the variable unset the root is `./docking_data`, relative to the working directory,
and `.gitignore` blocks it. `python paths.py` prints where the current campaign is reading
and writing, which is the first thing to run when the dashboard says there is no result.
