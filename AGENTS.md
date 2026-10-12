# Build layout

- Keep each project's build outputs under that project's root, normally `build/<preset>/`.
- Do not build FEM_mfem into the shared dependency directory `E:\Dev\fem`.
- Build spdlog, nlohmann/json and muparser through project-local FetchContent by default. Their sources and build artifacts belong in `build/<preset>/_deps/`.
- The shared dependency directory contains installed SDKs/libraries and their metadata, not application binaries, validation results or project build helpers.
- Keep Windows build configuration and the build/run entry point in this project. Preserve the WSL source project and existing simulation outputs when validating Windows builds.
- Manage Python tools with uv and the project's `pyproject.toml`/`uv.lock`. Keep the environment in `.venv/`, managed Python in `.uv/python/`, and its cache in `build/uv-cache/`. Do not put project Python environments in the shared native SDK directory.
