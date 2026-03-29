#include "fem/solver/LinearSolverFactory.hpp"

#include "fem/solver/MfemPcgSolver.hpp"
#include "fem/solver/MfemUmfpackSolver.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace fem::solver
{
namespace
{
std::string NormalizeSolverName(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name;
}
}  // namespace

std::unique_ptr<ILinearSolver> CreateLinearSolver(const std::string &solver_name, double rel_tol,
                                                  double abs_tol, int max_iter, int print_level)
{
    const std::string normalized = NormalizeSolverName(solver_name);
    if (normalized.empty() || normalized == "pcg")
    {
        return std::make_unique<MfemPcgSolver>(rel_tol, abs_tol, max_iter, print_level);
    }
    if (normalized == "umfpack")
    {
        return std::make_unique<MfemUmfpackSolver>();
    }

    throw std::runtime_error("Unsupported solver type: '" + solver_name +
                             "'. Supported solvers: pcg, umfpack.");
}
}  // namespace fem::solver
