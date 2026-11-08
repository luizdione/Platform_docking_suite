@echo off
title Docking Suite Dashboard
REM =============================================================================
REM Inicia o dashboard da plataforma no ambiente Python corrente.
REM
REM Nao ha caminho de maquina aqui. O script roda a partir da propria pasta e
REM usa o `python` do PATH. Se voce usa conda, ative o ambiente ANTES de chamar
REM este .bat, ou aponte DOCKING_PYTHON para o interprete desejado:
REM
REM     set DOCKING_PYTHON=C:\caminho\para\envs\meuenv\python.exe
REM     start_docking.bat
REM
REM Onde ficam os resultados: defina DOCKING_DATA_DIR. Sem ela, a plataforma
REM escreve em .\docking_data (relativo a esta pasta).
REM =============================================================================

cd /d "%~dp0"

if defined DOCKING_PYTHON (
    set "PY=%DOCKING_PYTHON%"
) else (
    set "PY=python"
)

"%PY%" --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo [ERRO] Python nao encontrado.
    echo Ative o seu ambiente, ou defina DOCKING_PYTHON com o caminho do interprete.
    pause
    exit /b 1
)

echo.
"%PY%" --version
if defined DOCKING_DATA_DIR (
    echo Raiz de dados: %DOCKING_DATA_DIR%
) else (
    echo Raiz de dados: .\docking_data   ^(defina DOCKING_DATA_DIR para mudar^)
)
echo.

"%PY%" dashboard_v2.py
