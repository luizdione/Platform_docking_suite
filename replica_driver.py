#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replica_driver.py - Roda REPLICAS INDEPENDENTES de MD sequencialmente, robusto a
quedas, por dias sem supervisao.

Por que replicas
----------------
Uma trajetoria unica (n=1) nao sustenta afirmacao sobre estabilidade de pose: o
que se mede nela pode ser o acaso de um sorteio de velocidades. Cada replica
reusa o MESMO pipeline (md_simulation.py), mas com TAG NOVO -> as velocidades
sao reatribuidas aleatoriamente no warm-up (setVelocitiesToTemperature sem
seed) e o integrador Langevin Middle e estocastico => replicas independentes a
partir da MESMA pose inicial.

Robustez
--------
- Subprocesso por corrida: um segfault do OpenMM mata so o filho; o driver segue.
- RESUME automatico: se uma corrida cai com .chk + minimized_<tag>.pdb gravados,
  a retentativa entra em modo --fast --complex minimized_<tag>.pdb, e md_simulation
  retoma a producao do checkpoint (mesmo nº de atomos -> casa). Se nao ha checkpoint
  ainda, retenta --fresh do zero.
- Estado idempotente em replica_state.json: re-lancar o driver pula o que ja terminou.
- Log mestre (replica_driver.log) + log por replica (md_<tag>_run.log).

Configuracao (nada de caminho de maquina soldado)
-------------------------------------------------
    DOCKING_DATA_DIR   raiz dos dados (ver paths.py)
    DOCKING_PYTHON     interprete das sub-rodadas (default: o que roda este script)
    DOCKING_MD_NS      nanossegundos por replica (default 50)
    DOCKING_MD_PLATFORM plataforma OpenMM (default OpenCL)

As replicas vem da linha de comando ou de um JSON:
    python replica_driver.py --target TcTR --replicas 3
    python replica_driver.py --jobs minhas_replicas.json
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys
import time
import datetime
from pathlib import Path

# ----------------------------------------------------------------------------
CODE  = Path(__file__).resolve().parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
import paths as _paths                                           # noqa: E402

PY    = os.environ.get("DOCKING_PYTHON") or sys.executable
MDSIM = CODE / "md_simulation.py"
DATA  = _paths.data_root()
WORK  = DATA / "figures" / "md_clean"
STATE = DATA / "replica_state.json"
STATUS = DATA / "replica_status.json"
MASTER = DATA / "replica_driver.log"
AGG   = DATA / "replica_summary.json"

NS = float(os.environ.get("DOCKING_MD_NS", "50"))
PLATFORM = os.environ.get("DOCKING_MD_PLATFORM", "OpenCL")
MAX_RETRIES = 2          # tentativas extras por corrida (resume)

# Preenchido por main() a partir dos argumentos: nenhuma lista de corridas de um
# estudo especifico fica soldada aqui.
JOBS = []

# Ambiente dos subprocessos. Se o interprete vem de um env conda, o prefixo e
# deduzido dele -- antes era um caminho absoluto de uma maquina so.
ENV = os.environ.copy()
_envdir = Path(PY).resolve().parent
if (_envdir / "Library" / "bin").exists():                       # layout conda/Windows
    ENV["CONDA_PREFIX"] = str(_envdir)
    ENV["PATH"] = os.pathsep.join([
        str(_envdir),
        str(_envdir / "Library" / "bin"),
        str(_envdir / "Library" / "mingw-w64" / "bin"),
        str(_envdir / "Scripts"),
        ENV.get("PATH", ""),
    ])
ENV.setdefault("OPENMM_DEFAULT_PLATFORM", PLATFORM)


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    line = "[%s] %s" % (now(), msg)
    print(line, flush=True)
    try:
        with open(MASTER, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:
        pass


def load_json(p, default):
    try:
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        pass
    return default


def save_json(p, obj):
    try:
        p.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    except Exception as e:
        log("  [WARN] nao salvou %s: %s" % (p.name, e))


def md_done(tag):
    """True se a producao terminou: resumo RMSD existe com t_total ~50 ns."""
    j = WORK / ("md_%s_rmsd_summary.json" % tag)
    if not j.exists():
        return False
    try:
        d = json.loads(j.read_text(encoding="utf-8"))
        return float(d.get("t_total_ns", 0)) >= 49.0
    except Exception:
        return False


def build_cmd(tag, target):
    """Escolhe --fresh (1a vez) ou --fast+resume (se ha checkpoint do mesmo sistema)."""
    chk = WORK / ("md_%s.chk" % tag)
    minpdb = WORK / ("minimized_%s.pdb" % tag)
    base = [PY, str(MDSIM),
            "--target", target, "--mode", "subversive",
            "--ns", str(NS), "--tag", tag,
            "--platform", PLATFORM, "--allow-fallback"]
    if chk.exists() and minpdb.exists():
        # retoma do checkpoint usando exatamente o sistema solvatado minimizado
        # (mesmo nº de atomos -> loadCheckpoint casa). md_simulation faz resume.
        return base + ["--fast", "--complex", str(minpdb)], "resume(--fast)"
    return base + ["--fresh"], "fresh"


def run_md_job(tag, target):
    cmd, mode = build_cmd(tag, target)
    runlog = WORK / ("md_%s_run.log" % tag)
    log("  [%s] modo=%s :: %s" % (tag, mode, " ".join(cmd)))
    t0 = time.time()
    try:
        with open(runlog, "a", encoding="utf-8") as fh:
            fh.write("\n\n==== %s :: %s ====\n" % (now(), mode))
            fh.flush()
            rc = subprocess.run(cmd, cwd=str(DOCK), env=ENV,
                                stdout=fh, stderr=subprocess.STDOUT).returncode
    except Exception as e:
        log("  [%s] EXCECAO no subprocess: %s" % (tag, e))
        rc = -99
    dt = (time.time() - t0) / 3600.0
    return rc, dt


def _build_jobs(args):
    """Monta a lista de replicas: --jobs <json> ou --target/--replicas."""
    if args.jobs:
        data = json.loads(pathlib.Path(args.jobs).read_text(encoding="utf-8"))
        jobs = data["jobs"] if isinstance(data, dict) else data
        for j in jobs:
            if not j.get("target") or not j.get("tag"):
                raise SystemExit("cada corrida precisa de 'target' e 'tag': %r" % j)
        return jobs
    if not args.target:
        raise SystemExit("informe --target (com --replicas) ou --jobs <arquivo.json>")
    first = args.first_replica
    return [{"target": args.target, "tag": "%s_%s_r%d" % (args.target, args.mode, n)}
            for n in range(first, first + args.replicas)]


def main():
    ap = argparse.ArgumentParser(
        description="Roda replicas independentes de MD em serie, com resume.")
    ap.add_argument("--target", help="id do alvo no registry (ex.: TcTR)")
    ap.add_argument("--replicas", type=int, default=3, help="quantas replicas (default 3)")
    ap.add_argument("--first-replica", type=int, default=1, dest="first_replica",
                    help="numero da primeira replica, para continuar uma serie")
    ap.add_argument("--mode", default="run", help="rotulo do modo, usado no tag")
    ap.add_argument("--jobs", help="JSON com [{target, tag}, ...] para controle total")
    args = ap.parse_args()

    global JOBS
    JOBS = _build_jobs(args)

    log("=" * 72)
    log("replica_driver INICIADO | %d corridas | plataforma=%s | %g ns cada"
        % (len(JOBS), PLATFORM, NS))
    log("python=%s" % PY)
    log("workdir=%s" % WORK)
    log("ordem: " + " -> ".join(j["tag"] for j in JOBS))
    log("=" * 72)

    if not MDSIM.exists():
        log("FATAL: md_simulation.py nao encontrado em %s" % MDSIM)
        sys.exit(2)

    state = load_json(STATE, {})
    t_start = time.time()

    for i, job in enumerate(JOBS, 1):
        tag, target = job["tag"], job["target"]
        write_status(state, tag, t_start)

        if state.get(tag, {}).get("status") == "done" or md_done(tag):
            log("[%d/%d] %s: ja concluido; pulando." % (i, len(JOBS), tag))
            state[tag] = {"status": "done"}
            save_json(STATE, state)
            continue

        attempts = int(state.get(tag, {}).get("attempts", 0))
        ok = False
        while attempts <= MAX_RETRIES:
            attempts += 1
            log("[%d/%d] %s: tentativa %d/%d" % (i, len(JOBS), tag, attempts, MAX_RETRIES + 1))
            state[tag] = {"status": "running", "attempts": attempts, "started": now()}
            save_json(STATE, state)
            write_status(state, tag, t_start)

            rc, dt = run_md_job(tag, target)
            done = md_done(tag)
            log("[%d/%d] %s: rc=%s  md_done=%s  (%.1f h)" % (i, len(JOBS), tag, rc, done, dt))

            if rc == 0 and done:
                log("[%d/%d] %s: MD OK." % (i, len(JOBS), tag))
                state[tag] = {"status": "done", "hours": round(dt, 2),
                              "attempts": attempts, "finished": now()}
                save_json(STATE, state)
                ok = True
                break
            else:
                last = attempts > MAX_RETRIES
                state[tag] = {"status": ("failed" if last else "retry"),
                              "attempts": attempts, "rc": rc, "last": now()}
                save_json(STATE, state)
                log("[%d/%d] %s: %s" % (i, len(JOBS), tag,
                    "DESISTINDO deste tag (segue p/ o proximo)" if last
                    else "vai retentar (resume se houver checkpoint)"))

        write_status(state, tag, t_start)
        elapsed = (time.time() - t_start) / 3600.0
        log("  >>> %s: %s | decorrido total: %.1f h"
            % (tag, "OK" if ok else "FALHOU", elapsed))
        # agrega parcialmente a cada corrida (deixa um overview sempre fresco)
        try:
            aggregate()
        except Exception as e:
            log("  [WARN] aggregate parcial falhou: %s" % e)

    log("=" * 72)
    log("TODAS AS CORRIDAS PROCESSADAS. Resumo final:")
    for j in JOBS:
        t = j["tag"]
        log("   %-22s %s" % (t, state.get(t, {}).get("status", "?")))
    aggregate()
    write_status(state, None, t_start)
    log("FIM replica_driver. Total: %.1f h" % ((time.time() - t_start) / 3600.0))


if __name__ == "__main__":
    main()
