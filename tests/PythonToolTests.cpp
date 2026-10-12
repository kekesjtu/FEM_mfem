#include "app/PythonTool.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;
using fem::app::PythonCommand;
using fem::app::ResolvePython;
using fem::app::RunPython;
using fem::app::SplitToolArguments;

void Check(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void SetEnvironment(const char *name, const std::string &value)
{
#ifdef _WIN32
    Check(_wputenv_s(fs::u8path(name).c_str(), fs::u8path(value).c_str()) == 0, "set environment");
#else
    Check(setenv(name, value.c_str(), 1) == 0, "set environment");
#endif
}

fs::path VenvPython(const fs::path &root)
{
#ifdef _WIN32
    return root / "Scripts/python.exe";
#else
    return root / "bin/python";
#endif
}

void Marker(const fs::path &path)
{
    fs::create_directories(path.parent_path());
    std::ofstream(path) << "resolver test marker\n";
#ifndef _WIN32
    fs::permissions(path, fs::perms::owner_exec, fs::perm_options::add);
#endif
}

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
            throw std::runtime_error("Usage: fem_python_tool_tests PYTHON_EXECUTABLE ARTIFACT_DIR");
        const PythonCommand actual{fs::absolute(fs::u8path(argv[1])), {}};
        const auto root = fs::absolute(fs::u8path(argv[2]));
        fs::create_directories(root);

        const std::vector<std::string> parsed{
            "--plot-dir", "folder with spaces", "--fields", "V,T", "--empty", "",
            "--path", "C:\\plain\\file.txt", "--text", "a \"quote\""};
        Check(SplitToolArguments(R"args(--plot-dir "folder with spaces" --fields 'V,T' --empty "" --path C:\plain\file.txt --text "a \"quote\"")args") == parsed,
              "quoted arguments and Windows path separators");
        bool rejected = false;
        try { SplitToolArguments("--plot-dir \"unclosed"); }
        catch (const std::runtime_error &) { rejected = true; }
        Check(rejected, "unclosed quotes must be rejected");

        const auto local = VenvPython(root / "project/.venv");
        const auto active = VenvPython(root / "active environment");
        Marker(local);
        Marker(active);
        SetEnvironment("FEM_PYTHON_EXECUTABLE", "");
        SetEnvironment("VIRTUAL_ENV", (root / "active environment").u8string());
        Check(ResolvePython(root / "project").executable == local, "project venv takes priority");
        Check(ResolvePython(root / "without-local").executable == active, "activated venv fallback");
        SetEnvironment("FEM_PYTHON_EXECUTABLE", actual.executable.u8string());
        Check(ResolvePython(root / "project").executable == actual.executable, "explicit override wins");
        SetEnvironment("FEM_PYTHON_EXECUTABLE", (root / "missing-python").u8string());
        rejected = false;
        try { ResolvePython(root / "project"); }
        catch (const std::runtime_error &) { rejected = true; }
        Check(rejected, "invalid explicit override must not silently fall back");

        SetEnvironment("FEM_PYTHON_EXECUTABLE", "");
        SetEnvironment("VIRTUAL_ENV", "");
#ifdef _WIN32
        const auto path_python = root / "path with spaces/python.exe";
#else
        const auto path_python = root / "path with spaces/python3";
#endif
        Marker(path_python);
        SetEnvironment("PATH", path_python.parent_path().u8string());
        Check(ResolvePython(root / "without-local").executable == path_python, "system PATH fallback");
#ifdef _WIN32
        Marker(root / "launcher/py.exe");
        SetEnvironment("PATH", (root / "launcher").u8string());
        const auto launcher = ResolvePython(root / "without-local");
        Check(launcher.executable == root / "launcher/py.exe" &&
                  launcher.prefix_arguments == std::vector<std::string>{"-3"}, "Windows Python launcher");
#endif
        SetEnvironment("PATH", "");
        rejected = false;
        try { ResolvePython(root / "without-local"); }
        catch (const std::runtime_error &) { rejected = true; }
        Check(rejected, "missing Python must produce an actionable error");

        const std::string expected = R"json(["two words", "雪", "\"quoted\"", "C:\\tail\\", "%PATH%", "$(echo x)", "&", ""])json";
        std::vector<std::string> launch_args{
            "-c", "import json,sys; assert sys.argv[2:] == json.loads(sys.argv[1]), repr(sys.argv); sys.exit(7)",
            expected, "two words", "雪", "\"quoted\"", "C:\\tail\\", "%PATH%", "$(echo x)", "&", ""};
        Check(RunPython(actual, launch_args) == 7, "literal argv and child exit code");
        const auto spaced_venv = root / fs::u8path("解释器 with spaces");
        Check(RunPython(actual, {"-m", "venv", "--without-pip", spaced_venv.u8string()}) == 0,
              "create interpreter at a Unicode path with spaces");
        Check(RunPython({VenvPython(spaced_venv), {}}, launch_args) == 7,
              "launch interpreter at a Unicode path with spaces");
        rejected = false;
        try { RunPython({root / "nonexistent", {}}, {"-c", "pass"}); }
        catch (const std::system_error &) { rejected = true; }
        Check(rejected, "process creation failure must be reported");
        std::cout << "PASS: Python discovery, quoted/literal arguments, Unicode paths, exit codes and errors\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
