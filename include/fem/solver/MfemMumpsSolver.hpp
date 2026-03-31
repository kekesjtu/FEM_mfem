#pragma once

#include "fem/solver/ILinearSolver.hpp"

namespace fem::solver
{
class MfemMumpsSolver final : public ILinearSolver
{
  public:
    MfemMumpsSolver() = default;

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;
};
}  // namespace fem::solver
