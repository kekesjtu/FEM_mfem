#pragma once

#include "fem/solver/ILinearSolver.hpp"

namespace fem::solver
{
class MfemAmgSolver final : public ILinearSolver
{
  public:
    MfemAmgSolver(double rel_tol, double abs_tol, int max_iter, int print_level);

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;

  private:
    double rel_tol_;
    double abs_tol_;
    int max_iter_;
    int print_level_;
};
}  // namespace fem::solver
