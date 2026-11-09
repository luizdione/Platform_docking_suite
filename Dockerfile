# =============================================================================
# Docking Suite — imagem de execucao da plataforma.
# A imagem carrega SO o codigo e as dependencias. Nenhum dado entra nela (ver
# .dockerignore). Monte a sua arvore de dados num volume e aponte
# DOCKING_DATA_DIR para o ponto de montagem.
#
# Base CUDA (nao conda-forge): o requirements.txt do projeto e 100% pip-instalavel
# (openmm>=8.4 publica wheel PyPI com suporte CUDA desde a serie 8.x; RDKit e
# openbabel-wheel idem). Isso replica o metodo "ALTERNATIVO (pip)" ja documentado
# no proprio requirements.txt — nao precisa de conda aqui.
# =============================================================================
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONUTF8=1 \
    PIP_NO_CACHE_DIR=1 \
    PATH="/opt/venv/bin:$PATH"

# --- Python 3.11 + libs de sistema exigidas por RDKit/OpenBabel/OpenMM ---------
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3-pip \
        wget unzip ca-certificates \
        libxrender1 libxext6 libsm6 build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

WORKDIR /app

# --- Criação do Ambiente Virtual Isolado para mitigar erros de dependência -----
RUN python3 -m venv /opt/venv

# --- Dependencias Python (camada cacheavel) ------------------------------------
COPY requirements.txt .

# Utiliza o resolver antigo emparelhado com o ambiente isolado do venv
RUN pip install --upgrade pip setuptools wheel && \
    pip install --use-deprecated=legacy-resolver -r requirements.txt

# --- AutoDock Vina (nao instalavel via pip — baixado do release oficial) ------
# Binario Linux oficial do AutoDock Vina (ccsb-scripps). O valor anterior baixava
# a HOME do github.com como "vina" -- a imagem nunca rodaria docking. --version
# no build FALHA a imagem se o binario nao for executavel (sanidade fail-loud).
ARG VINA_VERSION=1.2.5
RUN wget -qL "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
        -O /usr/local/bin/vina && \
    chmod +x /usr/local/bin/vina && \
    /usr/local/bin/vina --version

# --- Metadados OCI ------------------------------------------------------------
# `image.source` e o que faz o GitHub LIGAR o pacote a este repositorio
# automaticamente, e com isso herdar a visibilidade dele a cada push. Sem o
# label, o pacote nasce orfao e privado, e o vinculo tem de ser feito a mao na
# pagina do pacote -- foi o que aconteceu na primeira publicacao.
LABEL org.opencontainers.image.source="https://github.com/luizdione/Platform_docking_suite" \
      org.opencontainers.image.title="Docking Suite" \
      org.opencontainers.image.description="Umbrella platform for structure-based antiparasitic screening. Code and configuration only - no research results." \
      org.opencontainers.image.authors="Luiz Dione Barbosa De-Melo <luiz.melo@ifrj.edu.br>" \
      org.opencontainers.image.url="https://github.com/luizdione/Platform_docking_suite" \
      org.opencontainers.image.documentation="https://github.com/luizdione/Platform_docking_suite#readme" \
      org.opencontainers.image.licenses="LicenseRef-DockingSuite-1.0"


# --- Codigo da plataforma -------------------------------------------------------
COPY . .

# Raiz de dados DENTRO do container. /data e um ponto de montagem vazio na
# imagem: sem volume, a plataforma escreve ali e o que for escrito morre com o
# container -- de proposito, para que nada de dado persista na imagem.
#   docker run -e DOCKING_DATA_DIR=/data -v /host/dados:/data ...
ENV DOCKING_DATA_DIR=/data
RUN mkdir -p /data

# Porta do dashboard Flask (dashboard_v2.py aceita --host/--port)
EXPOSE 5001

# ENTRYPOINT dual:
#   - modo job: 1 execucao de cada vez via runner.py (CLI com
#     proveniencia: plan/run/status/resume).
#   - modo dashboard: inspecao manual e relatorio.
# Ex.: docker run --gpus all IMG python3 runner.py run --campaign DEMO01 --kind docking
ENTRYPOINT ["python3"]
CMD ["dashboard_v2.py", "--host", "0.0.0.0", "--port", "5001"]
