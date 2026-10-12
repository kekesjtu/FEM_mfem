#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace fem::app
{
struct PythonCommand
{
    std::filesystem::path executable;
    std::vector<std::string> prefix_arguments;
};

// FEM_PYTHON_EXECUTABLE, project .venv, VIRTUAL_ENV, then PATH.
PythonCommand ResolvePython(const std::filesystem::path &project_root);

// Single/double quotes and escaped quotes are supported; no shell expansion.
std::vector<std::string> SplitToolArguments(const std::string &arguments);

// Launch directly and return the child's exit status. Launch errors throw.
int RunPython(const PythonCommand &python, const std::vector<std::string> &arguments);
} // namespace fem::app
