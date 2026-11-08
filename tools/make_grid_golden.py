#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Grava o golden GEOMETRICO: o centro da caixa de docking de cada alvo.

Por que existe
--------------
`tests/test_registry_parity.py` e `tests/test_engines.py` comparam o centro de
grade calculado hoje contra um valor gravado antes. Esse e o gate que pega a
regressao silenciosa: trocar um `ligand_resnames`, um `grid_strategy` ou a
estrutura baixada move o centro da caixa, o docking continua "funcionando" e
devolve outro resultado sem erro nenhum.

O golden NAO vem no repositorio, de proposito: calcula-lo exige os receptores
preparados, que vivem na sua arvore de dados ($DOCKING_DATA_DIR), nao aqui. Sem
o golden os dois testes dizem SKIP com o motivo -- e e voce quem decide quando
congelar o valor.

Quando gravar
-------------
Depois de preparar os receptores de uma campanha e ANTES de rodar a campanha de
verdade. O que o golden congela e a geometria de entrada, nao o resultado.

Quando regravar
---------------
So quando a mudanca de centro for deliberada (voce trocou a estrutura de
referencia, o sitio ou a politica de caixa). Regravar para "fazer o teste
passar" devolve o gate ao estado de nao existir.

Uso
---
    python tools/make_grid_golden.py                 # campanha corrente
    python tools/make_grid_golden.py --campaign DEMO01
    python tools/make_grid_golden.py --dry-run       # mostra, nao grava
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import paths        # noqa: E402
import registry     # noqa: E402

GOLDEN = ROOT / "tests" / "golden" / "grid_centers_golden.json"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--campaign", default=None, help="id da campanha (default: a corrente)")
    ap.add_argument("--dry-run", action="store_true", help="imprime e nao grava")
    ap.add_argument("--out", default=None, help="caminho de saida (default: %s)" % GOLDEN.name)
    args = ap.parse_args()

    try:
        from Bio.PDB import PDBParser
    except ImportError:
        print("[ERRO] Biopython nao instalado (pip install biopython).", file=sys.stderr)
        return 2

    import docking_v2

    campaign = args.campaign or paths.default_campaign()
    protein_dir = paths.data_root(campaign) / "proteins"
    if not protein_dir.exists():
        print(f"[ERRO] receptores nao encontrados em {protein_dir}", file=sys.stderr)
        print("       Rode o docking da campanha uma vez, ou defina DOCKING_DATA_DIR.",
              file=sys.stderr)
        return 1

    db = registry.as_legacy_target_db(campaign)
    parser = PDBParser(QUIET=True)
    centers, missing = {}, []

    for tid, info in sorted(db.items()):
        pdb_id = info.get("pdb")
        pdb_path = protein_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            missing.append(f"{tid} ({pdb_id})")
            continue
        structure = parser.get_structure("prot", str(pdb_path))
        center = docking_v2._calculate_grid_center(structure, info)
        centers[tid] = {
            "pdb": pdb_id,
            "grid_center": [round(float(c), 2) for c in center],
            "grid_strategy": info.get("grid_strategy"),
            "grid_size": list(info.get("grid_size") or []),
        }
        print(f"  {tid:12s} {pdb_id:6s} {centers[tid]['grid_center']}")

    if missing:
        print(f"\n[AVISO] sem PDB para: {', '.join(missing)}", file=sys.stderr)
    if not centers:
        print("[ERRO] nenhum centro calculado -- nada a gravar.", file=sys.stderr)
        return 1

    payload = {
        "_note": ("Centro da caixa de docking por alvo, congelado a partir dos "
                  "receptores preparados. Gate contra regressao silenciosa de "
                  "geometria: se o centro mudar, o resultado de docking muda."),
        "_campaign": campaign,
        "_generated_by": "tools/make_grid_golden.py",
        "_protein_dir": str(protein_dir),
        "centers": centers,
    }
    out = Path(args.out) if args.out else GOLDEN
    if args.dry_run:
        print("\n[dry-run] nao gravado. Conteudo:")
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\n[OK] {len(centers)} centro(s) gravados em {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
