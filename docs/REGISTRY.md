# The declarative registry

The single source of truth about **what** is docked, **against what**, and **under which
protocol**. Code reads the registry; code contains no domain dictionary.

Before it existed, the platform had three target dictionaries that disagreed with each
other — one of them assigned to *TcCruzain* the PDB code `3GRS`, which is human
glutathione reductase. The registry exists so that cannot happen again.

## Hierarchy

```
registry/
  program.json                              umbrella: identity, chemical defaults
  campaigns/<ID>.json                       campaign (= one study, = one panel)
  targets/<id>.json                         TARGET   — global, reusable
  ligands/<id>.json                         LIGAND   — global, reusable
  pairs/<CAMPAIGN>__<target>__<ligand>.json PAIR     — scoped to a campaign
  protocols/<id>.json                       PROTOCOL — docking, structure, MD
```

**Targets and ligands are global; the pair belongs to the campaign.** The same hGR is a
selectivity counter-target in one study and could be the primary target in another: the
target does not change, the pair does.

It was the absence of the pair concept that let `md_simulation.py` reuse
`cofactor_FAD_TcTR.xml` for targets that are not flavoenzymes.

## The pair is where the peculiarity lives

Everything that exists only *in the combination* of target × ligand, and that used to sit
implicit in a script, is a required field here:

| Field | Why it belongs to the pair, not to the target or the ligand |
|---|---|
| `assembly_state.required` | The jasplakinolide site exists only at the junction of three protomers of the filament; the latrunculin site exists in the monomer. Same target, opposite requirements. |
| `experimental_kd[].state` | **Required.** Cytochalasin D is 4.6 µM on G-actin and ~2 nM at the barbed end. A Kd without its state is a category error. |
| `engines.affinity_estimator_reliable` | An estimator can be trustworthy on one system and useless on another. That is configuration, not an assumption. |
| `validation.training_leakage` | If the reference PDB predates the predictor's training cutoff, a hit proves nothing. |
| `validation.verdict_rule` | Reading only `model_0` of a co-folding run is an order statistic dressed up as an answer. |
| `cofactor_requirement` | FAD in trypanothione reductase and ATP/Mg in actin are not accessories: without them the site ceases to exist. |

## Fields that change pipeline behaviour

**Target** — `structure.source` (`pdb` \| `alphafold` \| `local`). A predicted model
requires `confidence.global` (pLDDT), and the registry refuses the file without it.
`chemistry.cofactor_source` (`crystal` \| `transferred` \| `none`); `transferred` requires
`cofactor_transfer.{method,donor_pdb,align_scope}` — the superposition used to transfer a
cofactor may need to be **local to the pocket**, not global, when the target has
interdomain motion.

**Site** — `sites.<name>.definition` accepts five modes. Three were already implemented by
`_calculate_grid_center()` (`ligand_center`, `metal_center`, `residue_center`); two were
added to cover a target with no co-crystal of its own
(`derived_from_reference_complex`, `explicit_center`). `box_policy` says whether the box is
fixed or derived from the bounding box of the crystallographic ligand.

**Ligand** — `topology.is_macrocycle` forces `prep.pdbqt.rigid_macrocycles: false`: a
19-membered ring treated as rigid produces a wrong pose by construction, with no error.
`prep.protonation.required` forces `post_check`, because protonation fails silently.

## Use

```python
import registry

registry.get_target("TcTR")                  # accepts id or alias ("1BZL", "3GRS")
registry.resolve_site("TcTR")                # default site, with its id resolved
registry.get_pair("DEMO01", "TcTR", "beta-lapachone")
registry.resolve_protocol("DEMO01", "md", "TcTR", "beta-lapachone")   # pair > campaign
registry.list_targets("DEMO01")
registry.validate_all()                      # [] = everything is fine
```

`DOCKING_REGISTRY` overrides the directory; `DOCKING_CAMPAIGN` overrides the dashboard's
current campaign.

## Writing rules

- **Explicit UTF-8, always.** `notes` fields carry scientific prose, and the Windows
  default (cp1252) corrupts it silently.
- **Generate, do not transcribe.** The targets and ligands came out of
  `tools/export_target_db_to_registry.py`. Transcribing 10 targets × 12 fields by hand is
  exactly how the wrong PDB codes were born.
- **English.** These files are public documentation, not internal notes.
- **A missing field stays missing.** Do not fill one in with a plausible value.

## Validation

```bash
python registry.py                      # summary + validation
python registry.py --campaign DEMO01    # the campaign's targets
python -m pytest tests/test_registry_parity.py
```

`validate_all()` is a hand-written validator, deliberately: `jsonschema` is not in
`requirements.txt`, and the Windows/conda environment is fragile enough already. If
`jsonschema` happens to be installed and `registry/_schema/` holds files, they are used as
an additional check.

## Bridge to legacy code

`registry.as_legacy_target_db(campaign)` projects the registry into the shape of the old
`docking_v2.TARGET_DB`, with the same types — including `grid_size` as a **tuple**, because
a list unpacks identically (`sx, sy, sz = grid_size`) and the regression would pass
unnoticed at runtime. `docking_v2.load_target_db()` consumes it and falls back, loudly, to
an embedded table if the registry disappears.
