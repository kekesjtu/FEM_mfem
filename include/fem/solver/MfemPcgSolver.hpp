#pragma once

#include "fem/solver/ILinearSolver.hpp"

#include <memory>

namespace fem::solver
{
class MfemPcgSolver final : public ILinearSolver
{
  public:
    MfemPcgSolver(double rel_tol, double abs_tol, int max_iter, int print_level);

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;
    void Mult(const mfem::Vector &b, mfem::Vector &x) override;

  private:
    double rel_tol_;
    double abs_tol_;
    int max_iter_;
    int print_level_;
    std::unique_ptr<mfem::Solver> preconditioner_;
    std::unique_ptr<mfem::CGSolver> solver_;
};
}  // namespace fem::solver
