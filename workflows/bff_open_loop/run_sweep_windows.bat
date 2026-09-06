@echo off
title X-56 BFF open-loop sweep
wsl.exe bash -lc "cd /home/nicomonzi/X_56/workflows/bff_open_loop && ./run_sweep.sh"
if errorlevel 1 (
  echo.
  echo SWEEP FALLITO. Controllare i file .stdout nella cartella dei risultati.
) else (
  echo.
  echo SWEEP COMPLETATO.
)
pause
