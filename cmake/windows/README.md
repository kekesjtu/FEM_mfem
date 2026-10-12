# Native Windows build

Build from the project root:

```bat
tools\windows.cmd build
```

The command discovers the installed MSVC Build Tools, initializes its x64
environment, and uses the `windows-msvc-intel-full-release` CMake preset.
Its default SDK directory is `E:\Dev\fem`; set `FEM_WINDOWS_SDK_ROOT` to use
another installed SDK directory.

Project outputs belong in this project's own tree:

```text
build/windows-msvc-intel-full-release/
  bin/FEM_solution.exe
  lib/
  mfem/
  _deps/spdlog-src, spdlog-build
  _deps/nlohmann_json-src, nlohmann_json-build
  _deps/muparser-src, muparser-build
```

spdlog 1.15.3, nlohmann/json 3.11.3 and muparser 2.3.5 are obtained and built
by FetchContent. `FEM_USE_SYSTEM_LIGHTWEIGHT_DEPENDENCIES=ON` is an explicit
opt-in to installed packages; it is OFF by default.

The shared SDK directory contains installed hypre/METIS/MUMPS/MFEM/OpenMP
libraries and Intel MPI/oneMKL/ifx. It does not contain application
builds or project build/run scripts.

Python comparison tools use the uv-managed environment in the project root.
The dependencies are declared in `pyproject.toml` and locked in `uv.lock`.
Restore them with:

```bat
tools\windows.cmd python-sync
```

This keeps `.venv/`, the managed Python installation (`.uv/python/`), and
the download cache (`build/uv-cache/`) inside this project. The `python3`
command used by the C++ comparison hook resolves to this environment.
For standalone tool invocation, use `uv run tools/compare_solution_txt.py ...`
from the project root.

Run a configuration from the project root:

```bat
set OMP_NUM_THREADS=2
tools\windows.cmd run 2 path/to/config.json
```

Relative mesh and output paths in the configuration are interpreted relative
to the project root. Validation configurations redirect output into the build
directory so existing simulation results are preserved.

All C/C++ code uses MSVC. The Fortran portion of MUMPS uses ifx. The current
MFEM revision requires `/openmp:llvm` and UTF-8 source decoding; the installed
OpenMP runtime is selected by the run command.
