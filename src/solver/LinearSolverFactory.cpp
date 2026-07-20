#include "solver/LinearSolverFactory.hpp"

#include "solver/MfemPcgSolver.hpp"
#if defined(FEM_ENABLE_MUMPS)
#include "solver/MfemMumpsSolver.hpp"
#endif
#if defined(FEM_ENABLE_CPARDISO)
#include "solver/MfemPardisoSolver.hpp"
#endif

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
                                                  double abs_tol, int max_iter, int print_level,
                                                  mfem::ParFiniteElementSpace *elasticity_fespace)
{
    const std::string normalized = NormalizeSolverName(solver_name);
    if (normalized.empty() || normalized == "pcg")
    {
        return std::make_unique<MfemPcgSolver>(rel_tol, abs_tol, max_iter, print_level,
                                               elasticity_fespace);
    }
#if defined(FEM_ENABLE_MUMPS)
    if (normalized == "mumps")
    {
        return std::make_unique<MfemMumpsSolver>();
    }
#endif
#if defined(FEM_ENABLE_CPARDISO)
    if (normalized == "pardiso")
    {
        return std::make_unique<MfemPardisoSolver>();
    }
#endif

    throw std::runtime_error("Unsupported solver type: '" + solver_name +
                             "'. Supported solvers: pcg"
#if defined(FEM_ENABLE_MUMPS)
                             ", mumps"
#endif
#if defined(FEM_ENABLE_CPARDISO)
                             ", pardiso"
#endif
                             ".");
}
}  // namespace fem::solver
