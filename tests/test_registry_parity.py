# -*- coding: utf-8 -*-
"""
Gate de paridade dos alvos.

Roda em TODAS as fases da generalizacao multi-alvo. Hoje compara a fonte de
verdade vigente (docking_v2.TARGET_DB) contra os goldens da Fase 0; a partir da
Fase 2, docking_v2.TARGET_DB passa a ser servido pelo registry declarativo e
este mesmo teste vira a prova de que a troca foi transparente.

Se algum destes testes falhar, a refatoracao mudou ciencia -- nao so codigo.

Uso:
    pytest tests/test_registry_parity.py -v
    python tests/test_registry_parity.py          # sem pytest instalado
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

# Receptores preparados: vivem na arvore de DADOS, nao no repositorio.
# Sem eles o gate geometrico nao tem o que medir e se declara skip -- o que e
# honesto: e um teste de integracao, nao de unidade.
def _protein_dir():
    import paths
    return paths.data_root() / "proteins"

# Campos cujo valor define o resultado cientifico. Mudanca aqui = mudanca de
# resultado. Campos de prosa (notes, literature) ficam de fora de proposito:
# podem ser reescritos sem invalidar nada.
SCIENTIFIC_FIELDS = [
    "pdb",
    "organism",
    "target",
    "category",
    "keep_cofactors",
    "ligand_resnames",
    "grid_strategy",
    "grid_fallback_residues",
    "grid_size",
    "human_ortholog",
]


def _load(name):
    with open(GOLDEN_DIR / name, encoding="utf-8") as fh:
        return json.load(fh)


def _current_target_db():
    import docking_v2

    return docking_v2.TARGET_DB


def test_target_set_unchanged():
    """Nenhum alvo do golden pode sumir nem trocar de nome."""
    golden = _load("demo01_targets_golden.json")["targets"]
    current = _current_target_db()
    assert set(current) >= set(golden), (
        "alvos sumiram da fonte de verdade: "
        f"{sorted(set(golden) - set(current))}"
    )


def test_scientific_fields_unchanged():
    """Todo campo que define resultado permanece identico ao golden."""
    golden = _load("demo01_targets_golden.json")["targets"]
    current = _current_target_db()
    diffs = []

    for tid, gold in golden.items():
        cur = current.get(tid)
        if cur is None:
            diffs.append(f"{tid}: ausente")
            continue
        for field in SCIENTIFIC_FIELDS:
            got, want = cur.get(field), gold.get(field)
            # tupla vs lista: a ida-e-volta pelo JSON converte grid_size.
            # Comparamos o CONTEUDO, e o tipo e checado separadamente.
            if isinstance(got, tuple):
                got = list(got)
            if isinstance(want, tuple):
                want = list(want)
            if got != want:
                diffs.append(f"{tid}.{field}: golden={want!r} atual={got!r}")

    assert not diffs, "campos cientificos divergiram do golden:\n  " + "\n  ".join(diffs)


def test_grid_size_is_tuple():
    """grid_size tem que continuar tupla depois de passar pelo registry JSON.

    Uma lista desempacota igual (sx, sy, sz = grid_size), entao a regressao
    passaria despercebida em runtime -- mas quebra comparacoes de igualdade.
    """
    current = _current_target_db()
    offenders = [
        tid for tid, info in current.items() if not isinstance(info["grid_size"], tuple)
    ]
    assert not offenders, (
        f"grid_size deixou de ser tupla em: {offenders} "
        "(provavel ida-e-volta pelo JSON sem re-tuplificar)"
    )


def test_grid_centers_unchanged():
    """Gate geometrico: o centro da caixa de docking de cada alvo.

    Exige os receptores preparados na arvore de dados e um golden gravado numa
    execucao anterior. Sem um dos dois, SKIP com o motivo -- um gate que
    "passa" sem ter medido nada e pior que um gate ausente.
    """
    import pytest
    from Bio.PDB import PDBParser

    import docking_v2

    golden_path = GOLDEN_DIR / "grid_centers_golden.json"
    if not golden_path.exists():
        pytest.skip(
            "sem tests/golden/grid_centers_golden.json. Gere-o depois de preparar "
            "os receptores: python tools/make_grid_golden.py")
    golden = json.loads(golden_path.read_text(encoding="utf-8"))["centers"]
    protein_dir = _protein_dir()
    if not protein_dir.exists():
        pytest.skip(f"receptores ausentes em {protein_dir} (defina DOCKING_DATA_DIR)")
    parser = PDBParser(QUIET=True)
    current = _current_target_db()
    diffs, skipped = [], []

    for tid, gold in golden.items():
        info = current.get(tid)
        if info is None:
            diffs.append(f"{tid}: alvo ausente")
            continue
        pdb_path = protein_dir / f"{gold['pdb']}.pdb"
        if not pdb_path.exists():
            skipped.append(tid)
            continue
        structure = parser.get_structure("prot", str(pdb_path))
        center = [round(float(c), 2) for c in docking_v2._calculate_grid_center(structure, info)]
        if center != gold["grid_center"]:
            diffs.append(f"{tid}: golden={gold['grid_center']} atual={center}")

    assert not diffs, "CENTRO DE GRADE MUDOU (invalida o docking):\n  " + "\n  ".join(diffs)
    if len(skipped) > 1:
        pytest.skip(f"PDBs demais ausentes para validar: {skipped}")


def test_no_duplicate_pdb_ids():
    """Guarda para o lookup reverso por PDB em docking_v2.clean_protein().

    `[k for k,v in TARGET_DB.items() if v['pdb']==...][0]` devolve o alvo errado
    se dois alvos compartilharem PDB, e estoura IndexError se o alvo nao tiver
    PDB (caso AlphaFold, que chega com a campanha de actina). Enquanto esse
    lookup existir, esta invariante precisa valer.
    """
    current = _current_target_db()
    seen = {}
    for tid, info in current.items():
        pdb = info.get("pdb")
        assert pdb, f"{tid} nao tem PDB -- o lookup reverso em clean_protein() vai estourar"
        seen.setdefault(pdb, []).append(tid)
    dups = {p: t for p, t in seen.items() if len(t) > 1}
    assert not dups, f"PDB compartilhado entre alvos, lookup reverso fica ambiguo: {dups}"


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"[PASS] {name}")
        except AssertionError as exc:
            failures += 1
            print(f"[FAIL] {name}\n       {exc}")
    print(f"\n{failures} falha(s)")
    sys.exit(1 if failures else 0)
