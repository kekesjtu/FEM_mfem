#pragma once

#include "fem/solver/ILinearSolver.hpp"

namespace fem::solver
{
class MfemUmfpackSolver final : public ILinearSolver
{
  public:
    MfemUmfpackSolver() = default;

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;
};
}  // namespace fem::solver
