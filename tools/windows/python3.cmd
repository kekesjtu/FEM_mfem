@echo off
setlocal EnableExtensions
for %%I in ("%~dp0..\..") do set "FEM_PROJECT_ROOT=%%~fI"
if not exist "%FEM_PROJECT_ROOT%\.venv\Scripts\python.exe" (
  echo Project Python environment missing. Run tools\windows.cmd python-sync.
  exit /b 1
)
"%FEM_PROJECT_ROOT%\.venv\Scripts\python.exe" %*
exit /b %errorlevel%
