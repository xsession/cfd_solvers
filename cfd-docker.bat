@echo off
REM cfd_solvers Docker helper launcher.
REM All logic lives in cfd-docker.ps1 (same directory) - no bash dependency.
REM
REM Usage:
REM   cfd-docker.bat build
REM   cfd-docker.bat run [case] [--threads N]
REM   cfd-docker.bat test
REM   cfd-docker.bat mpi [--ranks N]
REM   cfd-docker.bat help
REM
REM Environment: CFD_TAG (image tag), CFD_MPI=1 (build MPI image).

setlocal
set "SCRIPT_DIR=%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%SCRIPT_DIR%cfd-docker.ps1" %*
set "EXIT=%ERRORLEVEL%"
endlocal & exit /b %EXIT%
