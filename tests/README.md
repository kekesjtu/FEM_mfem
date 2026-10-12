# Python tool regression checks

`PythonToolTests.cpp` covers interpreter selection (override, project/active
venv, PATH, and the Windows launcher), quoted CLI arguments, literal shell
metacharacters, Unicode/space paths, child exit codes, and launch failures.
It uses only the standard library and does not link to MFEM or MPI.

After configuring the main project, build the excluded test target in the same
compiler/SDK environment used for the project:

```text
cmake --build build/windows-msvc-intel-full-release --target fem_python_tool_tests
build/windows-msvc-intel-full-release/bin/fem_python_tool_tests.exe .venv/Scripts/python.exe build/python-tool-tests/windows-artifacts
```

The same sources can be checked separately on Linux/macOS:

```sh
mkdir -p build/python-tool-tests
c++ -std=c++17 -Wall -Wextra -Iinclude src/app/PythonTool.cpp tests/PythonToolTests.cpp -o build/python-tool-tests/linux-tests
build/python-tool-tests/linux-tests /usr/bin/python3 build/python-tool-tests/linux-artifacts
```

Tests create files and a temporary venv inside the supplied artifact directory.
That directory should be under the project's `build/` directory. Environment
changes affect only the test process.
