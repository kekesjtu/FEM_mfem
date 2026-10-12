@echo off
setlocal EnableExtensions
for %%I in ("%~dp0..") do set "FEM_PROJECT_ROOT=%%~fI"
if not defined FEM_WINDOWS_SDK_ROOT set "FEM_WINDOWS_SDK_ROOT=E:\Dev\fem"
set "FEM_ONEAPI=%FEM_WINDOWS_SDK_ROOT%\oneapi\Library"
set "PATH=%FEM_PROJECT_ROOT%\tools\windows;%FEM_PROJECT_ROOT%\.venv\Scripts;%FEM_WINDOWS_SDK_ROOT%\install\openmp\bin;%FEM_ONEAPI%\bin;%FEM_ONEAPI%\bin\libfabric\utils;%PATH%"
set "I_MPI_ROOT=%FEM_ONEAPI%"
set "MKLROOT=%FEM_ONEAPI%"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
set "MPLBACKEND=Agg"
set "MPLCONFIGDIR=%FEM_PROJECT_ROOT%\build\matplotlib-cache"
set "UV_PROJECT_ENVIRONMENT=%FEM_PROJECT_ROOT%\.venv"
set "UV_PYTHON_INSTALL_DIR=%FEM_PROJECT_ROOT%\.uv\python"
set "UV_CACHE_DIR=%FEM_PROJECT_ROOT%\build\uv-cache"
set "UV_PYTHON_INSTALL_BIN=0"
set "UV_PYTHON_NO_REGISTRY=1"
cd /d "%FEM_PROJECT_ROOT%"
if /i "%~1"=="build" goto build
if /i "%~1"=="run" goto run
if /i "%~1"=="python-sync" goto python_sync
goto usage
:build
set "FEM_VS_DIR="
for /f "usebackq tokens=*" %%I in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "FEM_VS_DIR=%%I"
if not defined FEM_VS_DIR (
  echo MSVC Build Tools were not found.
  exit /b 1
)
call "%FEM_VS_DIR%\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64
if errorlevel 1 exit /b 1
set "PATH=%FEM_VS_DIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;%FEM_VS_DIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja;%PATH%"
set "LIB=%FEM_ONEAPI%\lib;%LIB%"
set "INCLUDE=%FEM_WINDOWS_SDK_ROOT%\oneapi\opt\compiler\include\intel64;%FEM_ONEAPI%\include;%INCLUDE%"
cmake --preset windows-msvc-intel-full-release
if errorlevel 1 exit /b 1
cmake --build --preset build-windows-msvc-intel-full-release --parallel 8
exit /b %errorlevel%
:run
if "%~2"=="" goto usage
if "%~3"=="" goto usage
set "FEM_EXE=%FEM_PROJECT_ROOT%\build\windows-msvc-intel-full-release\bin\FEM_solution.exe"
if not exist "%FEM_EXE%" (
  echo Run tools\windows.cmd build first.
  exit /b 1
)
mpiexec -localonly -n %2 "%FEM_EXE%" "%~3"
exit /b %errorlevel%
:python_sync
uv sync --locked --managed-python
exit /b %errorlevel%
:usage
echo Usage: tools\windows.cmd build
echo        tools\windows.cmd run MPI_RANKS CONFIG_JSON
echo        tools\windows.cmd python-sync
exit /b 2
