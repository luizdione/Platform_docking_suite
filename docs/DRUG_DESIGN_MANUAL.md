# Computational drug design module
## Docking Suite — IFRJ

**Module:** `drug_design_module.py`
**Dashboard tab:** Drug Design
**Parent scaffold:** whatever you give it. The default is β-lapachone, a published natural
product — no compound is hardcoded in the module.
**Central reference for the catalogue:** Pinto & De Castro, *Molecules* 2009;
De Castro et al., *Front. Physiol.* 2014

> **What this manual documents, and what it is not.** The R-group catalogue, the synthetic
> routes and the ΔG model all come from **published SAR** — each entry carries its citation
> in the `ref` field. Nothing here is a research result: the module proposes synthesis
> hypotheses from the scaffold you hand it, and the priority is always to validate them in
> real docking.

---

## What the module does

It generates **new antiparasitic candidates** in silico from the compounds in the
database, using three complementary strategies:

1. **FBDD** — systematic substitution of functional groups (R-groups) on the quinone scaffold
2. **Hybridisation** — fusion of BRICS fragments from two active compounds (the LQB-118 design strategy)
3. **Genetic algorithm** — evolution guided by a fitness function combining ADMET, estimated affinity and synthetic accessibility

Each generated candidate gets:
- a canonical SMILES validated by RDKit
- a computational ADMET profile (Lipinski, QED, TPSA, SAS)
- an estimated binding ΔG against TcTR (empirical QSAR model)
- a **suggested synthetic route**, based on established protocols
- a **scientific rationale** grounded in the cited literature
- a results paragraph ready for a report (EN/PT)

---

## Chemical basis — why the quinone scaffold

### The target: trypanothione reductase

Trypanosomatids do not use glutathione as their principal antioxidant thiol: they use
**trypanothione**, kept reduced by trypanothione reductase (TR). The pathway does not exist
in mammals, and the closest human counterpart is glutathione reductase (hGR), whose active
site is narrower. That difference in architecture is what leaves room for a selective
ligand — and it is also why the platform docks the parasite enzyme **always** alongside the
human counterpart. Without the counter-target, "binds well" is not selectivity.

### Why a quinone

Two properties, documented in the antiparasitic literature, make the 1,4-quinone core a
recurring starting point:

- **Redox cycling.** The quinone is reduced to a semiquinone, which reacts with O₂ and
  returns the quinone while generating superoxide. The result is sustained oxidative
  stress — and it is relatively selective for the parasite, whose thiol-dependent defence
  is TR itself.
- **Subversive substrate.** Quinones compete for electrons at the flavin face of
  flavoenzymes, diverting the flux that should have reduced the physiological thiol.

Hence the catalogue below: each R-group is a modification at C2 of the core with
**published** antiparasitic SAR, and the `ref` field says where it comes from.

> The two mechanisms are neither mutually exclusive nor automatic. Which one dominates
> depends on the molecule and the target, and that is an experimental question — the module
> generates the hypothesis, not the answer.

---

## Approach 1 — FBDD (fragment-based drug design)

### Principle

It replaces the functional group at **C2 of 1,4-naphthoquinone** (between the two
carbonyls) with R-groups drawn from a 20-group catalogue based on published SAR.

```
         O                     O
         ‖                     ‖
    ┌────C────┐           ┌────C────┐
    │         │           │         │
  ──C         C─── [X]    C         C─── [new R-group]
    │         │           │         │
    └────C────┘           └────C────┘
         ‖                     ‖
         O                     O

    Parent compound          FBDD analogue
```

### R-group catalogue, by category

#### Halogens — redox cycling and halogen bonds

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–Cl** | Cl₂/AcOH or NCS/DMF | ↑logP +0.7; halogen bond to the backbone; NQ1-Cl active | Salomão et al. 2013 BMC Microbiol |
| **–Br** | Br₂/AcOH or NBS/EtOAc | Greater polarisability than Cl; NQ1-Br IC50 ~0.8 μM vs T. cruzi | Duarte et al. 2020 |
| **–F** | SELECTFLUOR® or F-TEDA-BF₄ | Bioisostere of H; C–F resists CYP450; ↑metabolic stability | Rani 2022 Chem Biol Drug Des |

#### Sulfur groups — interaction with the Cys of the TcTR site

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–SOCH₃** | mCPBA (1 equiv., 0 °C, DCM) | ↑solubility vs –SCH₃; sulfoxides are metabolically stable | Rani 2022 review |
| **–S-thiazolyl** | SNAr + 2-mercaptothiazole/K₂CO₃/DMF | π-stacking with the FAD of TcTR; imidazole bioisostere | thiazolyl-NQ SAR |
| **–SO₂NH₂** | Oxone®/MeOH:H₂O | TPSA +42 Å²; H-bond to Glu18 and Tyr110 of TcTR | De Castro group reviews |

#### Nitrogen groups — H-bonds and activation by TcNTR

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–NH₂** | SnCl₂/EtOH/HCl or NH₃/autoclave | H-bond to Glu18; the 2-amino-NQ series is well studied | Pinto & De Castro 2009 Molecules |
| **–NO₂** | conc. HNO₃/H₂SO₄ (0 °C) | Dual action: ROS + activation by TcNTR, as in benznidazole | Docampo & De Souza 1978; Salomão 2013 |
| **–triazole** | CuAAC click chemistry | IC50 0.8–3.1 μM, more active than benznidazole; easy synthesis, >90% yield | Bahia et al. 2016 MedChemComm |
| **–imidazole** | SNAr of 2-halo-NQ + imidazole/K₂CO₃/DMF | Basis of the N1/N2/N3 naphthoimidazoles, IC50 0.5–2.0 μM | Pinto & De Castro 2009 Molecules |
| **–aminopyrazine** | Buchwald-Hartwig or SNAr with Cl-pyrazine | Basic N → ↑solubility at acidic pH (the Leishmania endosome) | aminopyrazine medicinal chemistry |

#### Oxygen groups — solubility and bioisosteres

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–OH** (lawsone-like) | SeO₂/dioxane or K₂CO₃/MeOH hydrolysis | LQB-118 precursor; Fe-SOD chelation; TPSA +17 Å² | Maia 2011 Invest New Drugs |
| **–OCH₃** | Me₂SO₄/K₂CO₃/acetone | Blocks oxidation; logP +0.5; atovaquone analogue | Lopes 1978 Ann Trop Med |
| **–OCF₃** | AgOCF₃/DCM or Umemoto | ↑↑lipophilicity; CYP450 resistance; used in modern antiprotozoals | OCF₃ medicinal chemistry |

#### Fused heterocycles — π-stacking with FAD

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–isoxazolyl** | [3+2] cycloaddition, metal-free | Imidazole bioisostere; modern antiparasitics | MDPI Pharmaceutics 2022 |
| **–pyrazolyl** | Hydrazine + β-ketoester condensation | pyrazolyl-NQ analogues active vs T. brucei | Rani 2022 review |

#### Aliphatic groups — lapachol/menadione-like

| R-group | Synthesis | Rationale | Reference |
|---|---|---|---|
| **–CH₃** (menadione) | MeMgBr/CuI/THF | Vitamin K3; active vs T. cruzi | Docampo et al. 1978 Z Parasitenkd |
| **–CF₃** | Togni I/II or TMSCF₃/CsF | σp +0.54; resists metabolism; used in modern antiprotozoals | CF₃ SAR in antiparasitics |
| **–allyl** (lapachol) | Allylamine/K₂CO₃/DMF | Lapachol fragment; Michael acceptor in vivo | Moura 2001 J Braz Chem Soc |
| **–dihydrobenzofuranyl** | Oxa-Michael + cyclisation (the LQB-118 synthesis) | Inspired by LQB-118 (a pterocarpanquinone) | Maia 2011; Cunha-Junior 2016 AAC |

### How to use it in the dashboard

1. Click **FBDD (R-groups)**
2. Enter the SMILES of the parent compound, or click **Load active compound**
3. Select R-group categories (Ctrl+click for several)
4. Set the maximum number of variants (20 by default)
5. Click **Generate candidates**
6. Read the ranked candidate table
7. Click **View synthesis** on any candidate for the full detail

### Reading the FBDD results

| Metric | Good | Moderate | Poor |
|---|---|---|---|
| Estimated ΔG | ≤ −8.0 kcal/mol | −6.0 to −8.0 | > −6.0 |
| Drug score | ≥ 7.0 | 5.0 – 7.0 | < 5.0 |
| SAS | ≤ 3.0 (easy) | 3.0 – 5.0 | > 5.0 |
| Lipinski | pass | — | fail |
| TPSA | 40–90 Å² | 90–140 Å² | > 140 Å² |

---

## Approach 2 — molecular hybridisation (BRICS)

### Principle — after LQB-118

**Molecular hybridisation** is a design strategy widely used in the antiparasitic
literature. **LQB-118** (a pterocarpanquinone) was created by fusing **lapachol** (a
natural naphthoquinone) with a **pterocarpan** (a plant fragment), giving:

- IC50 = 2.5 μM against T. cruzi epimastigotes
- 99% reduction of parasite load in vivo (visceral leishmaniasis, 10 mg/kg oral)
- activity in chronic myeloid leukaemia — an anticancer effect

```
     Lapachol                 Pterocarpan               LQB-118
  [NQ scaffold]    +    [bicyclic fragment]  →  [pterocarpanquinone hybrid]
  IC50 ~20 μM            — antifungal              IC50 2.5 μM T. cruzi
                                                   SI > 40x (CC50 >100 μM VERO)
```

### The module's BRICS algorithm

1. **Fragmentation**: each compound is broken at the bonds defined by the BRICS algorithm (Breaking Retrosynthetically Interesting Chemical Substructures — Degen et al. 2008)
2. **Filtering**: fragments with fewer than 4 heavy atoms are discarded
3. **Recombination**: every fragment-A + fragment-B combination is tried
4. **Lipinski filter**: MW ≤ 650, HBD ≤ 5
5. **Scoring and ranking**: drug_score + estimated ΔG + SAS

### Automatically suggested synthetic routes

The module reads the chemistry of the fragments and suggests the most suitable route:

| Situation | Suggested route |
|---|---|
| Fragment B contains N | Buchwald-Hartwig (Pd₂(dba)₃/Xphos) or SNAr |
| Fragment B contains O | Williamson condensation or oxa-Michael |
| Both aromatic | Suzuki-Miyaura (Pd/boronate) or Ullmann (Cu) |
| No heteroatoms | Knoevenagel or Michael addition |

### Example combinations

Pairs formed only from **published** compounds, among those that ship with the platform.
They let you see the module working before you use your own scaffold.

| Compound A | Compound B | Rationale |
|---|---|---|
| β-Lapachone | LQB-118 | Furanone ring (ROS + NQO1 bioactivation) + pterocarpan fragment |
| Plumbagin | β-Lapachone | Chelating –OH (active vs *Leishmania* Fe-SOD) + NQO1 bioactivation |
| β-Lapachone | Benznidazole | ROS + nitroreduction hybrid (two approved-drug mechanisms) |
| Lapachol | Atovaquone | Natural core + the lipophilic chain of a complex III ligand |

### How to use it in the dashboard

1. Click **Hybridisation**
2. Enter the SMILES of compound A (the parent)
3. Enter the SMILES of compound B (the partner, from the compound database)
4. Set the maximum number of hybrids
5. Click **Generate candidates**

---

## Approach 3 — genetic algorithm

### Darwinian evolution applied to chemistry

```
Generation 0:      Generation 5:      Generation 15:
[compound bank] → [mutations +]   → [optimised]
                   [crossovers]      [candidates]

Fitness = f(ΔG_TcTR) × f(ADMET) × f(SAS)
         ↑ affinity    ↑ drug-like  ↑ easy to make
```

### Genetic operations implemented

**Mutation (35% probability):** picks an R-group category at random and applies the
substitution to the parent compound via FBDD.

**Crossover (40% probability):** fragments two compounds via BRICS and recombines
fragments from different parents (hybridisation).

**Elitism (20% of the population):** the best individuals of each generation are carried
over unchanged.

**Tournament selection:** 3 random individuals compete; the fittest goes to reproduction.

### Fitness function

```
fitness = dg_score × admet_score × sas_score × lipinski_bonus

where:
  dg_score    = max(0, (-ΔG_est - 2.0) / 8.0)     [0–1, ΔG from -2 to -10 kcal/mol]
  admet_score = drug_score / 10.0                  [0–1]
  sas_score   = 1 - (SAS - 1.0) / 9.0              [0–1, SAS from 1 to 10]
  lip_bonus   = 1.0 if Lipinski passes, 0.4 if it fails
```

### Recommended parameters by goal

| Goal | Generations | Population | Estimated time |
|---|---|---|---|
| Quick exploration / smoke test | 10 | 15 | ~1–2 min |
| Standard screen | 15 | 20 | ~3–5 min |
| Detailed optimisation | 25 | 30 | ~10–15 min |
| Design campaign | 50 | 50 | ~30–60 min |

### How to use it in the dashboard

1. Click **Genetic algorithm**
2. Enter the SMILES of the seed compounds, one per line — three from the bank is a good start
3. Enter the matching names, in the same order
4. Set generations and population size
5. Click **Generate candidates**
6. Follow convergence on the generated evolution plot

---

## QSAR model — estimated ΔG for TcTR

The module uses an empirical QSAR model to estimate the binding affinity (ΔG) of
candidates against **T. cruzi trypanothione reductase (TcTR, PDB 1BZL)** without running
docking for every variant.

### Model calibration

The model was calibrated on ~50 naphthoquinones with published docking data:

| Publication | Compounds | ΔG range |
|---|---|---|
| Pinto & De Castro 2009 | 16 NQ derivatives | -4.5 to -9.2 kcal/mol |
| Pinto & De Castro 2009 | 12 lapachol derivatives | -5.0 to -8.8 kcal/mol |
| Salomão et al. 2013 BMC Microbiol | 8 NQ analogues | -5.5 to -8.5 kcal/mol |
| Rani et al. 2022 Chem Biol Drug Des | 14 natural NQs | -4.8 to -9.1 kcal/mol |

### The model equation

```
ΔG_est = -5.5                        (base: plain 1,4-NQ)
       + f(logP)                      (sweet spot 2.0–3.5: -1.2)
       + f(TPSA)                      (ideal 50–90 Å²: -0.8)
       + f(HBD) × (-0.35)            (per H-bond donor)
       + f(HBA) × (-0.18)            (up to 5 acceptors)
       + f(MW)                        (ideal 200–380 Da: -0.6)
       + N_arom × (-0.5)             (aromatic N: π-stacking with FAD)
       + S_presence × (-0.7)         (sulfur at C2: Cys52/57)
       + N_hal × (-0.35)             (halogens: halogen bonds)
       + NO2_presence × (-0.5)       (TcNTR activation)
```

### Limits of the QSAR model

The estimated ΔG is a **relative ranking**, not an absolute prediction. Expect errors of
±1.5 to 2.5 kcal/mol against real Vina docking. For the top five candidates it is
**strongly recommended** to run real docking from the dashboard's Pipeline tab.

---

## Synthetic accessibility score (SAS)

SAS (Ertl & Schuffenhauer 2009) measures how easy a molecule is to synthesise, on a scale
of 1 to 10.

| SAS | Reading | Example structure |
|---|---|---|
| 1.0 – 2.5 | Very easy — commercial reagents, 1–2 steps | Menadione, plumbagin |
| 2.5 – 3.5 | Easy — standard synthesis, 2–4 steps | Lapachol, juglone |
| 3.5 – 5.0 | Moderate — 4–6 steps, protection/deprotection | β-lapachone, LQB-118 |
| 5.0 – 7.0 | Complex — needs specialist expertise | N1/N2/N3 naphthoimidazoles |
| 7.0 – 10.0 | Very complex — total synthesis | Advanced polycyclic structures |

**The module's criterion:** `synth_feasible = True` when SAS ≤ 4.5, that is, synthesis
reachable in a conventional organic chemistry laboratory.

---

## Recommended workflow

### From hypothesis to prioritised candidate

```
Step 1: FBDD on your scaffold
   → Categories: nitrogen (triazole, imidazole) + sulfur + halogen
   → Max variants: 20
   → Goal: find 3-5 analogues with estimated ΔG <= -8.0 and SAS <= 4.0

Step 2: Hybridisation
   → scaffold × beta-lapachone   (ROS + NQO1 mechanism)
   → scaffold × LQB-118          (pterocarpan fragment)

Step 3: GA from the best of the bank
   → Seeds: your scaffold + beta-lapachone + plumbagin
   → 15 generations, pop=20
   → Goal: reach a non-obvious scaffold

Step 4: VALIDATE -- the step you do not skip
   → Top-3 from each approach → real docking in the Pipeline tab
   → Compare against a reference drug (e.g. benznidazole)
   → The QSAR ΔG is triage; the docking ΔG is what goes in a report

Step 5: Report
   → "Download" button in the Drug Design tab (EN/PT)

Step 6: Synthesis proposal
   → Prioritise SAS <= 3.5 and ΔG <= -8.0
   → Hand over to a partner synthesis group
```

> Step 4 exists because the QSAR ΔG is a triage estimate, calibrated on published data for
> one chemical class. Treating it as an affinity is the error this module most easily
> invites.

### Priority criteria for synthesis

A candidate should satisfy **at least 3 of 4** criteria to be prioritised:

- [ ] Estimated ΔG ≤ −7.5 kcal/mol (better than the predicted reference drug)
- [ ] Drug score ≥ 7.0 / 10
- [ ] SAS ≤ 4.0 (accessible synthesis)
- [ ] Lipinski pass (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)

---

## Chaining drug design with the other modules

The module is designed to work in sequence with the rest of the pipeline:

```
Drug Design (new SMILES)
       ↓
Compounds (add to the bank via "Add new compound")
       ↓
Pipeline (run Vina docking to validate the real ΔG)
       ↓
ADMET (full pharmacokinetic profile)
       ↓
MD Simulation (confirm complex stability)
       ↓
Toxicity (compute SI = CC50/IC50 from your own data)
       ↓
Report (integrated text from the modules)
```

To add a drug-design candidate to the bank:

1. Copy the candidate's SMILES (from the table column or the "View synthesis" panel)
2. Go to the **Compounds** tab → SMILES field → paste → fill in a name → **Add**
3. The compound appears in the selector in the top banner
4. Select it as the active compound and run any tab as usual

---

## Glossary

| Term | Definition |
|---|---|
| **FBDD** | Fragment-Based Drug Design |
| **BRICS** | Breaking Retrosynthetically Interesting Chemical Substructures — fragmentation algorithm (Degen 2008) |
| **SAR** | Structure–Activity Relationship |
| **ΔG** | Gibbs free energy of binding (kcal/mol) — the more negative, the higher the affinity |
| **ADMET** | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| **QED** | Quantitative Estimate of Drug-likeness (Bickerton 2012) — 0 to 1, higher is better |
| **SAS** | Synthetic Accessibility Score (Ertl 2009) — 1 to 10, lower is easier to make |
| **TPSA** | Topological Polar Surface Area (Å²); 40–90 Å² is the usual target for antiparasitics |
| **TcTR** | *T. cruzi* trypanothione reductase — the target of the example campaign (PDB 1BZL) |
| **TcNTR** | *T. cruzi* nitroreductase — activates nitro prodrugs (benznidazole, –NO₂ analogues) |
| **SI** | Selectivity Index = CC50(host) / IC50(parasite) |
| **CuAAC** | Copper-catalysed Azide–Alkyne Cycloaddition — "click chemistry" |
| **SNAr** | Nucleophilic Aromatic Substitution |
| **GA** | Genetic Algorithm |
| **Fitness** | The function combining ΔG_TcTR + ADMET + SAS that steers the GA |

---

## Principal references

1. De Castro SL et al. *The double-edged sword in pathogenic trypanosomatids: the pivotal role of mitochondria in oxidative stress and bioenergetics.* Front Physiol 2014.

2. Pinto AV, De Castro SL. *The trypanocidal activity of naphthoquinones: a review.* Molecules 2009;14:4570–4590.

3. Menna-Barreto RFS, Henriques-Pons A, Pinto AV et al. *Effect of a β-lapachone-derived naphthoimidazole on Trypanosoma cruzi.* J Antimicrob Chemother 2005;56:1034–1041.

4. Bahia MT, Barbosa CB, Jardim et al. *Molecular hybridization as a powerful tool toward multitarget quinoidal systems.* MedChemComm 2016;7:1555.

5. Salomão K, De Santana NA, Molina MT, De Castro SL, Menna-Barreto RFS. *Mitochondrial swelling and membrane potential collapse as primary evidence of the mode of action of naphthoquinones.* BMC Microbiol 2013;13:196.

6. Maia RC et al. *LQB-118, a pterocarpanquinone structurally related to lapachol: a novel class of agent with high apoptotic effect in chronic myeloid leukemia cells.* Invest New Drugs 2011;29:1143.

7. Rani I, Bhatt A. *Natural naphthoquinones and their derivatives as potential drug molecules against trypanosome parasites.* Chem Biol Drug Des 2022;100:699–724.

8. Degen J et al. *On the art of compiling and using "drug-like" chemical fragment spaces.* ChemMedChem 2008;3:1503.

9. Ertl P, Schuffenhauer A. *Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions.* J Cheminform 2009;1:8.
