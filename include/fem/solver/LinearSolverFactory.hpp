#pragma once

#include "fem/solver/ILinearSolver.hpp"

#include <memory>
#include <string>

namespace fem::solver
{
std::unique_ptr<ILinearSolver> CreateLinearSolver(const std::string &solver_name, double rel_tol,
                                                  double abs_tol, int max_iter, int print_level);
}  // namespace fem::solver
