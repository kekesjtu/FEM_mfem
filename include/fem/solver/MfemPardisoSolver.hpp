#pragma once

#include "fem/solver/ILinearSolver.hpp"

#include <memory>

namespace fem::solver
{
class MfemPardisoSolver final : public ILinearSolver
{
  public:
    MfemPardisoSolver() = default;

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;
    void Mult(const mfem::Vector &b, mfem::Vector &x) override;

  private:
    std::unique_ptr<mfem::Solver> solver_;
};
}  // namespace fem::solver
