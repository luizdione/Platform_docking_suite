#!/usr/bin/env python3
"""
=============================================================================
COMPOUND DATABASE MODULE
Manages the ligand library in docking_project/ligand/
Each compound lives in its own subdirectory:
  docking_project/ligand/{compound_name}/
    metadata.json   — name, smiles, formula, mw, heavy_atoms, ...
    {name}.pdbqt    — prepared ligand for Vina
    {name}.mol2     — 3D structure
    {name}.sdf      — 3D structure (optional)
    {name}.svg      — 2D structure image (optional)
=============================================================================
"""
import json
import math
import sys
from pathlib import Path
from datetime import datetime
import logging

log = logging.getLogger(__name__)

import os as _os
import sys as _sys

# Raiz de dados: resolvida por paths.py ($DOCKING_DATA_DIR, default ./docking_data).
# Nunca uma pasta dentro do repositorio -- resultado nao mora no codigo.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in _sys.path:
    _sys.path.insert(0, str(_HERE))
try:
    import paths as _paths
    BASE_DIR = _paths.data_root()
except Exception:                                      # registry ausente
    BASE_DIR = Path(_os.environ.get("DOCKING_DATA_DIR") or (Path.cwd() / "docking_data"))
LIGAND_DIR = BASE_DIR / "ligand"

# ─── Biblioteca de referencia ────────────────────────────────────────────────
# Nao existe mais lista de compostos soldada aqui. A fonte unica de verdade sobre
# ligantes e `registry/ligands/*.json` (ver docs/REGISTRY.md): um dicionario de
# dominio dentro do codigo era a terceira copia da mesma informacao, e as tres
# divergiam. `reference_library()` projeta o registry na forma que este modulo
# consome, e `populate_library()` materializa essas entradas no banco em disco.


def reference_library(campaign_id=None) -> list[dict]:
    """Ligantes de referencia declarados no registry, na forma do compound_db.

    Com `campaign_id`, so os ligantes daquela campanha; sem, todos os do registry.
    Devolve [] se o registry nao estiver disponivel -- o chamador decide se isso
    e fatal (fail-loud no dashboard) ou apenas vazio (CLI).
    """
    try:
        import registry
    except Exception as exc:                                # pragma: no cover
        log.warning("reference_library: registry indisponivel: %s", exc)
        return []
    try:
        entries = (registry.list_ligands(campaign_id) if campaign_id
                   else [registry.get_ligand(i) for i in registry._ids("ligands")])
    except Exception as exc:
        log.warning("reference_library: falha ao ler ligantes: %s", exc)
        return []
    out = []
    for e in entries:
        out.append({
            "name": e.get("name") or e.get("id"),
            "smiles": e.get("smiles"),
            "display_name": e.get("display_name") or e.get("name") or e.get("id"),
            "formula": e.get("formula") or "",
            "mw": e.get("mw_g_mol"),
            "heavy_atoms": e.get("heavy_atoms"),
            "charge": e.get("charge", 0),
            "notes": e.get("notes") or "",
            "registry_id": e.get("id"),
        })
    return [e for e in out if e.get("smiles")]


def default_compound(campaign_id=None) -> dict | None:
    """Metadados do ligante default da campanha, lidos do registry."""
    try:
        import registry
        campaign = registry.get_campaign(campaign_id or registry.load_program()
                                         .get("default_campaign"))
        lig_id = campaign.get("default_ligand")
        if not lig_id:
            return None
        for e in reference_library():
            if e.get("registry_id") == lig_id:
                return e
    except Exception as exc:
        log.warning("default_compound: %s", exc)
    return None


# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def _ligand_root(base_dir=None, campaign_id=None) -> Path:
    """Raiz da biblioteca de ligantes.

    BASE_DIR/LIGAND_DIR sao resolvidos no import, o que basta enquanto ha uma
    campanha so. Com duas campanhas servidas pelo mesmo processo (raizes de
    dados diferentes), um valor de modulo devolveria a arvore errada -- dai o
    parametro. Sem argumento, o comportamento e exatamente o de antes.
    """
    if base_dir:
        return Path(base_dir) / "ligand"
    if campaign_id:
        try:
            import paths
            return paths.data_root(campaign_id) / "ligand"
        except Exception:
            pass
    return LIGAND_DIR


def get_ligand_dir(compound_name: str, base_dir=None, campaign_id=None) -> Path:
    """Return the directory for a given compound."""
    return _ligand_root(base_dir, campaign_id) / _safe_name(compound_name)


def _safe_name(name) -> str:
    """Nome de diretorio seguro para o filesystem.

    Rede de defesa contra entradas invalidas vindas do frontend: objetos JS
    serializados chegavam como '[object Object]' e viravam a pasta espuria
    '_object_Object_'. Rejeita nao-strings, vazios e sentinelas de serializacao,
    devolvendo 'compound' como fallback neutro. Nomes validos (ex.: 'beta-lapachone')
    passam inalterados.
    """
    import re
    if not isinstance(name, str):
        return "compound"
    s = name.strip()
    if not s or s.lower() in ("[object object]", "undefined", "null", "nan", "none"):
        return "compound"
    safe = re.sub(r'[^\w\-.]', '_', s)
    # rejeita resultados degenerados (so '_'/'.' ou o classico _object_Object_)
    if not safe.strip("_.") or safe.lower().startswith("_object_object"):
        return "compound"
    return safe


def list_compounds(base_dir=None, campaign_id=None) -> list[dict]:
    """
    Scan docking_project/ligand/ and return all compounds with their metadata.
    Returns list of metadata dicts, sorted by date_added then name.
    """
    compounds = []
    ligand_root = _ligand_root(base_dir, campaign_id)
    if not ligand_root.exists():
        return compounds
    for subdir in sorted(ligand_root.iterdir()):
        if not subdir.is_dir():
            continue
        meta_path = subdir / "metadata.json"
        if meta_path.exists():
            try:
                with open(meta_path, encoding='utf-8') as f:
                    meta = json.load(f)
                meta["_dir"] = str(subdir)
                meta["_has_pdbqt"] = any(subdir.glob("*.pdbqt"))
                meta["_has_mol2"]  = any(subdir.glob("*.mol2"))
                compounds.append(meta)
            except Exception as e:
                log.warning("list_compounds: metadata ilegivel em %s: %s", meta_path, e)
        else:
            # Directory exists but no metadata — infer minimal info
            compounds.append({
                "name":    subdir.name,
                "smiles":  "",
                "formula": "",
                "mw":      0,
                "_dir":    str(subdir),
                "_has_pdbqt": any(subdir.glob("*.pdbqt")),
                "_has_mol2":  any(subdir.glob("*.mol2")),
            })
    return compounds


def load_compound(name: str, base_dir=None, campaign_id=None) -> dict | None:
    """Load metadata for a specific compound. Returns None if not found."""
    cdir = get_ligand_dir(name, base_dir, campaign_id)
    meta_path = cdir / "metadata.json"
    if not meta_path.exists():
        # Try scanning all compounds by name field
        for c in list_compounds(base_dir, campaign_id):
            if c.get("name", "").lower() == name.lower():
                return c
        return None
    with open(meta_path, encoding='utf-8') as f:
        meta = json.load(f)
    meta["_dir"] = str(cdir)
    meta["_has_pdbqt"] = any(cdir.glob("*.pdbqt"))
    meta["_has_mol2"]  = any(cdir.glob("*.mol2"))
    return meta


def get_pdbqt_path(compound_name: str, base_dir=None, campaign_id=None) -> Path | None:
    """Return path to the .pdbqt file for a compound, or None if missing."""
    cdir = get_ligand_dir(compound_name, base_dir, campaign_id)
    pdbqt_files = list(cdir.glob("*.pdbqt"))
    if pdbqt_files:
        return pdbqt_files[0]
    return None


def add_compound(name: str, smiles: str, formula: str = "",
                 display_name: str = "", mw: float = 0.0,
                 notes: str = "", prepare_3d: bool = True) -> dict:
    """
    Add a new compound to the database.
    Creates the directory, saves metadata, and optionally generates 3D structure.
    Returns the metadata dict with status info.
    """
    safe = _safe_name(name)
    cdir = LIGAND_DIR / safe
    cdir.mkdir(parents=True, exist_ok=True)

    # Calculate MW from SMILES if not provided
    if mw == 0.0 or formula == "":
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if mw == 0.0:
                    mw = round(Descriptors.MolWt(mol), 3)
                if formula == "":
                    formula = rdMolDescriptors.CalcMolFormula(mol)
        except ImportError:
            log.debug("add_compound: RDKit ausente; MW/formula nao calculados de SMILES")

    meta = {
        "name":            name,
        "display_name":    display_name or name,
        "smiles":          smiles,
        "formula":         formula,
        "mw":              mw,
        "heavy_atoms":     _count_heavy_atoms(smiles),
        "rotatable_bonds": _count_rotatable_bonds(smiles),
        "charge":          0,
        "optimization":    "MMFF94 (1000 steps conjugate gradients)",
        "date_added":      datetime.now().strftime("%Y-%m-%d"),
        "notes":           notes,
    }

    # Save metadata
    with open(cdir / "metadata.json", "w", encoding='utf-8') as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    # Prepare 3D structure if requested
    status_msgs = []
    if prepare_3d:
        pdbqt_path = cdir / f"{safe}.pdbqt"
        mol2_path  = cdir / f"{safe}.mol2"
        if not pdbqt_path.exists():
            ok, msg = generate_pdbqt(smiles, pdbqt_path, mol2_path)
            status_msgs.append(msg)
        else:
            status_msgs.append(f"[OK] PDBQT already exists: {pdbqt_path.name}")

    # Generate 2D SVG structure image
    svg_path = cdir / f"{safe}.svg"
    if not svg_path.exists():
        _generate_svg(smiles, svg_path)

    meta["_dir"]        = str(cdir)
    meta["_has_pdbqt"]  = any(cdir.glob("*.pdbqt"))
    meta["_has_mol2"]   = any(cdir.glob("*.mol2"))
    meta["_status_msgs"] = status_msgs
    return meta


def generate_pdbqt(smiles: str, pdbqt_path: Path, mol2_path: Path = None) -> tuple[bool, str]:
    """
    Generate 3D structure (RDKit ETKDG+MMFF) and convert to PDBQT (obabel CLI).
    Falls back to pybel if obabel not found.  Returns (success, message).
    """
    # ── Método preferido: RDKit 3D + obabel CLI (robusto no env docking) ──
    try:
        import subprocess, tempfile
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("SMILES inválido")
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception as e:
            log.debug("generate_pdbqt: MMFF optimize falhou (%s); seguindo sem minimizar", e)
        pdb_tmp = str(pdbqt_path).replace(".pdbqt", "_3d.pdb")
        Chem.MolToPDBFile(mol, pdb_tmp)
        # obabel: PDB -> pdbqt
        obabel_exe = _find_obabel()
        r = subprocess.run([obabel_exe, pdb_tmp, "-opdbqt", "-O", str(pdbqt_path)],
                           capture_output=True, text=True, timeout=120)
        if Path(pdbqt_path).exists() and Path(pdbqt_path).stat().st_size > 0:
            if mol2_path:
                subprocess.run([obabel_exe, pdb_tmp, "-omol2", "-O", str(mol2_path)],
                               capture_output=True, text=True, timeout=60)
            return True, f"[OK] {pdbqt_path.name} (RDKit+obabel)"
        return False, f"[ERRO] obabel não gerou pdbqt: {r.stderr[:200]}"
    except Exception:
        pass  # fallback para pybel abaixo (RDKit/obabel indisponivel)

    # ── Fallback: pybel (pode falhar no conda se BABEL_LIBDIR não setado) ──
    try:
        from openbabel import pybel
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D(forcefield="mmff94", steps=500)

        # Energy minimisation
        ff = pybel._forcefields.get("mmff94")
        if ff and ff.Setup(mol.OBMol):
            ff.ConjugateGradients(1000, 1.0e-6)
            ff.GetCoordinates(mol.OBMol)

        mol.write("pdbqt", str(pdbqt_path), overwrite=True)
        if mol2_path:
            mol.write("mol2", str(mol2_path), overwrite=True)

        n_atoms = mol.OBMol.NumHvyAtoms()
        size    = pdbqt_path.stat().st_size if pdbqt_path.exists() else 0
        return True, f"[OK] 3D generated: {n_atoms} heavy atoms, {size:,} bytes"

    except ImportError:
        return False, "[ERRO] OpenBabel not installed (pip install openbabel-wheel)"
    except Exception as e:
        return False, f"[ERRO] 3D generation failed: {e}"


def _find_obabel() -> str:
    """Caminho do executavel obabel: $OBABEL_EXE, senao o do PATH."""
    import subprocess
    candidates = [c for c in (_os.environ.get("OBABEL_EXE"), "obabel") if c]
    for c in candidates:
        try:
            r = subprocess.run([c, "--version"], capture_output=True, timeout=5)
            if r.returncode == 0 or b"Open Babel" in (r.stdout + r.stderr):
                return c
        except Exception as e:
            log.debug("_find_obabel: candidato %r indisponivel: %s", c, e)
    return "obabel"  # fallback to PATH


def ensure_default_compound(campaign_id=None):
    """Garante que o ligante DEFAULT DA CAMPANHA exista no banco em disco.

    Qual composto e esse deixou de ser uma constante deste modulo: vem de
    `campaign.default_ligand` no registry. Assim a mesma plataforma serve
    campanhas com ligantes diferentes sem editar codigo.
    """
    entry = default_compound(campaign_id)
    if not entry:
        return False, ("[ERRO] campanha sem default_ligand no registry -- "
                       "declare registry/campaigns/<ID>.json:default_ligand")
    name = entry["name"]
    cdir = get_ligand_dir(name, campaign_id=campaign_id)
    meta_path = cdir / "metadata.json"

    if not meta_path.exists():
        cdir.mkdir(parents=True, exist_ok=True)
        with open(meta_path, "w", encoding='utf-8') as f:
            json.dump(entry, f, indent=2, ensure_ascii=False)

    safe = _safe_name(name)
    pdbqt = cdir / f"{safe}.pdbqt"
    if not pdbqt.exists():
        return generate_pdbqt(entry["smiles"], pdbqt, cdir / f"{safe}.mol2")
    return True, f"[OK] {name} ready"


def populate_library(overwrite: bool = False, campaign_id=None) -> list[dict]:
    """
    Materializa no banco em disco os ligantes de referencia do registry
    (`reference_library()`). Cria metadata.json para cada composto; nao
    sobrescreve entradas existentes a menos que overwrite=True. Nao gera 3D/PDBQT
    (prepare_3d=False) para ser rapido — rodar prep separadamente se necessario.
    Retorna lista de {name, status} para exibir no dashboard.
    """
    results = []
    library = reference_library(campaign_id)
    if not library:
        return [{"name": "-", "status": "registry sem ligantes -- nada a popular"}]
    for entry in library:
        name  = entry["name"]
        cdir  = get_ligand_dir(name)
        mpath = cdir / "metadata.json"
        if mpath.exists() and not overwrite:
            results.append({"name": name, "status": "já existe — pulado"})
            continue
        try:
            meta = add_compound(
                name=name,
                smiles=entry["smiles"],
                display_name=entry.get("display_name", name),
                notes=entry.get("notes", ""),
                prepare_3d=False,   # rápido; gerar 3D separadamente se necessário
            )
            results.append({"name": name, "status": "adicionado"})
        except Exception as e:
            results.append({"name": name, "status": f"ERRO: {e}"})
    return results


def clean_orphan_dirs() -> list[str]:
    """Remove subpastas em LIGAND_DIR sem metadata.json válido (ex: _object_Object_)."""
    removed = []
    if not LIGAND_DIR.exists():
        return removed
    for d in list(LIGAND_DIR.iterdir()):
        if not d.is_dir():
            continue
        mpath = d / "metadata.json"
        if not mpath.exists():
            try:
                import shutil
                shutil.rmtree(d)
                removed.append(d.name)
            except Exception as e:
                log.warning("clean_orphan_dirs: falha ao remover %s: %s", d, e)
    return removed


# =============================================================================
# SESSION MANAGEMENT
# =============================================================================

def make_session_tag() -> str:
    """Return a timestamp string for session subdirectories: YYYY-MM-DD_HH-MM-SS"""
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def get_session_dirs(base_dir: str | Path, session_tag: str, compound_name: str) -> dict:
    """
    Return a dict of timestamped output directories for a pipeline run.
    All per-session outputs live under figures/{compound}/{session}/ — a single
    canonical location.  The old docking/{compound}/{session}/ sub-tree was
    removed (it was never populated in practice and duplicated the layout).

    Structure:
      <base>/
        figures/          {compound}/{session}/   ← ALL per-session outputs
        reports/          {compound}/{session}/
        logs/             {compound}/{session}/
        results/          {compound}/             ← cumulative, not per-session
        proteins/                                 ← shared receptors
        ligand/           {compound}/             ← prepared ligand files
    """
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    sess   = session_tag

    dirs = {
        "ligand":   base / "ligand"   / safe_c,
        "proteins": base / "proteins",
        "results":  base / "results"  / safe_c,
        "figures":  base / "figures"  / safe_c / sess,
        "logs":     base / "logs"     / safe_c / sess,
        "texts":    base / "reports" / safe_c / sess,
        # "docking" alias kept for backward-compat callers; points to figures/
        "docking":  base / "figures"  / safe_c / sess,
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def list_sessions(base_dir: str | Path, compound_name: str) -> list[dict]:
    """
    List all past sessions for a compound, most recent first.
    Returns list of {tag, compound, figures_dir, docking_dir, texts_dir,
                     logs_dir, results_file}.
    Only looks in figures/ (single canonical location).
    """
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    sessions = []

    parent = base / "figures" / safe_c
    if parent.exists():
        for sess_dir in sorted(parent.iterdir(), reverse=True):
            if sess_dir.is_dir() and len(sess_dir.name) == 19:
                tag = sess_dir.name
                figures_dir = str(sess_dir)
                sessions.append({
                    "tag":          tag,
                    "compound":     compound_name,
                    "figures_dir":  figures_dir,
                    "docking_dir":  figures_dir,   # alias — same location
                    "texts_dir":    str(base / "reports" / safe_c / tag),
                    "logs_dir":     str(base / "logs" / safe_c / tag),
                    "results_file": str(base / "results" / safe_c / "results.json"),
                })
    return sorted(sessions, key=lambda x: x["tag"], reverse=True)


def get_latest_results(base_dir: str | Path, compound_name: str) -> dict | None:
    """Return the most recent results.json for a compound, or None."""
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    rf     = base / "results" / safe_c / "results.json"
    if rf.exists():
        # encoding explicito: sem ele o Windows abre em cp1252 e um resultado
        # com acento em qualquer campo de texto estoura UnicodeDecodeError
        with open(rf, encoding="utf-8") as f:
            return json.load(f)
    return None


# =============================================================================
# HELPERS
# =============================================================================

def _count_heavy_atoms(smiles: str) -> int:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol else 0
    except Exception:
        return 0


def _count_rotatable_bonds(smiles: str) -> int:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        return rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else 0
    except Exception:
        return 0


def _generate_svg(smiles: str, svg_path: Path) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem import Draw
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        mol = Draw.PrepareMolForDrawing(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_path.write_text(drawer.GetDrawingText())
        return True
    except Exception:
        return False


# =============================================================================
# CLI INTERFACE
# =============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compound Database Manager")
    sub = parser.add_subparsers(dest="cmd")

    # list
    sub.add_parser("list", help="List all compounds")

    # add
    p_add = sub.add_parser("add", help="Add a new compound")
    p_add.add_argument("--name",    required=True, help="Short ID (e.g. MY-COMPOUND-1)")
    p_add.add_argument("--smiles",  required=True)
    p_add.add_argument("--formula", default="")
    p_add.add_argument("--display", default="")
    p_add.add_argument("--notes",   default="")

    # info
    p_info = sub.add_parser("info", help="Show compound info")
    p_info.add_argument("name")

    # sessions
    p_sess = sub.add_parser("sessions", help="List past sessions for a compound")
    p_sess.add_argument("name")

    # init
    sub.add_parser("init", help="Garante o ligante default da campanha no banco")

    args = parser.parse_args()

    if args.cmd == "list":
        compounds = list_compounds()
        if not compounds:
            print("No compounds found. Run: python compound_db.py init")
        else:
            print(f"\n{'Name':<20} {'Formula':<16} {'MW':>8}  {'SMILES[:40]'}")
            print("-"*80)
            for c in compounds:
                smiles_preview = (c.get('smiles','') or '')[:40]
                pdbqt_ok = "[PDBQT]" if c.get("_has_pdbqt") else "[no PDBQT]"
                print(f"{c['name']:<20} {c.get('formula','?'):<16} {c.get('mw',0):>8.2f}  {smiles_preview} {pdbqt_ok}")

    elif args.cmd == "add":
        print(f"Adding compound: {args.name}")
        meta = add_compound(args.name, args.smiles, formula=args.formula,
                            display_name=args.display, notes=args.notes)
        print(f"  Directory: {meta['_dir']}")
        print(f"  Formula:   {meta['formula']}")
        print(f"  MW:        {meta['mw']}")
        print(f"  PDBQT:     {'OK' if meta['_has_pdbqt'] else 'MISSING'}")
        for msg in meta.get("_status_msgs", []):
            print(f"  {msg}")

    elif args.cmd == "info":
        c = load_compound(args.name)
        if c:
            for k, v in c.items():
                if not k.startswith("_"):
                    print(f"  {k:<20}: {v}")
        else:
            print(f"Compound '{args.name}' not found.")

    elif args.cmd == "sessions":
        sessions = list_sessions("docking_project", args.name)
        if not sessions:
            print(f"No sessions found for {args.name}")
        else:
            for s in sessions:
                print(f"  {s['tag']}  →  {s['figures_dir']}")

    elif args.cmd == "init":
        ok, msg = ensure_default_compound()
        print(msg)

    else:
        parser.print_help()
