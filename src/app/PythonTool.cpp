#include "app/PythonTool.hpp"

#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <system_error>

#ifdef _WIN32
#include <process.h>
#else
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif

namespace fem::app
{
namespace
{
namespace fs = std::filesystem;

std::string Environment(const char *name)
{
#ifdef _WIN32
    const auto wide_name = fs::u8path(name).wstring();
    wchar_t *value = nullptr;
    std::size_t count = 0;
    const auto error = _wdupenv_s(&value, &count, wide_name.c_str());
    if (error)
        throw std::system_error(error, std::generic_category(), "Could not read Python environment");
    const std::unique_ptr<wchar_t, decltype(&std::free)> storage(value, &std::free);
    return value ? fs::path(value).u8string() : std::string{};
#else
    const auto *value = std::getenv(name);
    return value ? value : std::string{};
#endif
}

bool IsExecutable(const fs::path &path)
{
    std::error_code error;
    if (!fs::is_regular_file(path, error))
        return false;
#ifdef _WIN32
    return true;
#else
    return access(path.c_str(), X_OK) == 0;
#endif
}

fs::path OnPath(const std::string &name)
{
    const auto path = fs::u8path(name);
    if (path.has_parent_path())
        return IsExecutable(path) ? fs::absolute(path) : fs::path{};
#ifdef _WIN32
    constexpr char separator = ';';
#else
    constexpr char separator = ':';
#endif
    const auto search = Environment("PATH");
    std::size_t start = 0;
    do
    {
        const auto end = search.find(separator, start);
        auto directory = search.substr(start, end - start);
        if (directory.size() >= 2 && directory.front() == '"' && directory.back() == '"')
            directory = directory.substr(1, directory.size() - 2);
        const auto candidate = fs::u8path(directory) / path;
        if (IsExecutable(candidate))
            return fs::absolute(candidate);
#ifdef _WIN32
        if (!path.has_extension() && IsExecutable(fs::path(candidate.wstring() + L".exe")))
            return fs::absolute(fs::path(candidate.wstring() + L".exe"));
#endif
        if (end == std::string::npos)
            break;
        start = end + 1;
    } while (start <= search.size());
    return {};
}

fs::path EnvironmentPython(const fs::path &root)
{
#ifdef _WIN32
    return root / "Scripts" / "python.exe";
#else
    return root / "bin" / "python";
#endif
}

#ifdef _WIN32
// CRT spawn joins argv with spaces; apply CRT quoting, including trailing '\\'.
std::wstring QuoteWindowsArgument(const std::string &argument)
{
    std::wstring result = L"\"";
    std::size_t slashes = 0;
    for (const auto ch : fs::u8path(argument).wstring())
    {
        if (ch == L'\\')
        {
            ++slashes;
            continue;
        }
        result.append(ch == L'"' ? slashes * 2 + 1 : slashes, L'\\');
        result += ch;
        slashes = 0;
    }
    result.append(slashes * 2, L'\\');
    result += L'"';
    return result;
}
#endif
} // namespace

PythonCommand ResolvePython(const fs::path &project_root)
{
    if (const auto configured = Environment("FEM_PYTHON_EXECUTABLE"); !configured.empty())
    {
        const auto executable = OnPath(configured);
        if (executable.empty())
            throw std::runtime_error("FEM_PYTHON_EXECUTABLE is not an executable: " + configured);
        return {executable, {}};
    }
    const auto local = EnvironmentPython(project_root / ".venv");
    if (IsExecutable(local))
        return {fs::absolute(local), {}};
    if (const auto active = Environment("VIRTUAL_ENV"); !active.empty())
    {
        const auto executable = EnvironmentPython(fs::u8path(active));
        if (IsExecutable(executable))
            return {fs::absolute(executable), {}};
    }
#ifdef _WIN32
    for (const auto *name : {"python.exe", "python3.exe"})
#else
    for (const auto *name : {"python3", "python"})
#endif
        if (const auto executable = OnPath(name); !executable.empty())
            return {executable, {}};
#ifdef _WIN32
    if (const auto launcher = OnPath("py.exe"); !launcher.empty())
        return {launcher, {"-3"}};
#endif
    throw std::runtime_error(
        "Python was not found. Restore the project's .venv with uv, activate a "
        "Python environment, or set FEM_PYTHON_EXECUTABLE to a Python 3 executable.");
}

std::vector<std::string> SplitToolArguments(const std::string &arguments)
{
    std::vector<std::string> result;
    std::string token;
    char quote = 0;
    bool started = false;
    for (std::size_t i = 0; i < arguments.size(); ++i)
    {
        const auto ch = arguments[i];
        if (ch == '\\' && quote != '\'' && i + 1 < arguments.size())
        {
            const auto next = arguments[i + 1];
            // Ordinary backslashes in Windows paths are preserved.
            if (next == '"' || next == '\\' ||
                (!quote && (next == '\'' || std::isspace(static_cast<unsigned char>(next)))))
            {
                token += next;
                ++i;
                started = true;
                continue;
            }
        }
        if (quote)
        {
            if (ch == quote)
                quote = 0;
            else
                token += ch;
        }
        else if (ch == '\'' || ch == '"')
        {
            quote = ch;
            started = true;
        }
        else if (std::isspace(static_cast<unsigned char>(ch)))
        {
            if (started)
            {
                result.push_back(token);
                token.clear();
                started = false;
            }
        }
        else
        {
            token += ch;
            started = true;
        }
    }
    if (quote)
        throw std::runtime_error("Unclosed quote in Python comparison arguments");
    if (started)
        result.push_back(token);
    return result;
}

int RunPython(const PythonCommand &python, const std::vector<std::string> &arguments)
{
    std::vector<std::string> values{python.executable.u8string()};
    values.insert(values.end(), python.prefix_arguments.begin(), python.prefix_arguments.end());
    values.insert(values.end(), arguments.begin(), arguments.end());
#ifdef _WIN32
    std::vector<std::wstring> quoted;
    for (const auto &value : values)
        quoted.push_back(QuoteWindowsArgument(value));
    std::vector<const wchar_t *> argv;
    for (const auto &value : quoted)
        argv.push_back(value.c_str());
    argv.push_back(nullptr);
    const auto status = _wspawnv(_P_WAIT, python.executable.c_str(), argv.data());
    if (status == -1)
        throw std::system_error(errno, std::generic_category(), "Could not launch Python");
    return static_cast<int>(status);
#else
    std::vector<char *> argv;
    for (auto &value : values)
        argv.push_back(value.data());
    argv.push_back(nullptr);
    pid_t child = 0;
    const auto error = posix_spawn(&child, python.executable.c_str(), nullptr, nullptr, argv.data(), environ);
    if (error)
        throw std::system_error(error, std::generic_category(), "Could not launch Python");
    int status = 0;
    while (waitpid(child, &status, 0) == -1)
        if (errno != EINTR)
            throw std::system_error(errno, std::generic_category(), "Could not wait for Python");
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
#endif
}
} // namespace fem::app
