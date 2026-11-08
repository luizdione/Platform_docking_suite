# -*- coding: utf-8 -*-
"""Gate de seguranca -- todo alvo OLIGOMERICO declara o estado de montagem e tem
centro de caixa EXPLICITO (nunca cai no centro de massa em silencio).

O sitio da colchicina na tubulina fica na interface intradimero: nao existe no
monomero, e o centro vem da POSE NATIVA do ligante no co-cristal
(definition=explicit_center, box_center_A gravado por
site_derivation.validate_embedded_ligand). Este teste garante que (a) os campos de
assembly estao vivos e validados, (b) todo sitio oligomerico tem box_center_A ->
docking_v2 nao tem como cair no ultimo recurso (centro de massa).

Uso:
    pytest tests/test_oligomer_site_safety.py -v
    python tests/test_oligomer_site_safety.py
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry  # noqa: E402

TARGETS = ROOT / "registry" / "targets"


def _t(tid):
    return json.loads((TARGETS / f"{tid}.json").read_text(encoding="utf-8"))


def test_tubulin_counter_target_colchicine_site():
    """Contra-alvo tubulina: categoria counter_target e sitio da colchicina explicito."""
    t = _t("tubulin_ab")
    assert t["category"] == "counter_target", t.get("category")
    s = t["structure"]
    assert (s.get("n_protomers_required") or 1) > 1 and s.get("assembly"), "assembly da tubulina"
    site = t["sites"]["colchicina"]
    assert site["definition"] == "explicit_center", site["definition"]
    assert isinstance(site.get("box_center_A"), list) and len(site["box_center_A"]) == 3


def test_every_multiprotomer_target_declares_assembly():
    """Invariante da validacao WS1: n_protomers_required>1 EXIGE assembly descrito."""
    for tid in registry._ids("targets"):
        s = (registry._read(registry._dir("targets") / f"{tid}.json").get("structure") or {})
        if (s.get("n_protomers_required") or 1) > 1:
            assert s.get("assembly"), f"{tid}: n_protomers>1 sem assembly (validacao WS1 falhou)"


def test_registry_valid_with_oligomers():
    """O registry inteiro valida sem erro."""
    errs = registry.validate_all()
    assert errs == [], "registry invalido:\n" + "\n".join(errs)


if __name__ == "__main__":
    failures = 0
    for _name, _fn in sorted(globals().items()):
        if not _name.startswith("test_") or not callable(_fn):
            continue
        try:
            _fn()
            print(f"[PASS] {_name}")
        except AssertionError as exc:
            failures += 1
            print(f"[FAIL] {_name}\n       {exc}")
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"[ERRO] {_name}\n       {type(exc).__name__}: {exc}")
    print(f"\n{failures} falha(s)")
    sys.exit(1 if failures else 0)
