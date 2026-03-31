#pragma once

#include "fem/solver/ILinearSolver.hpp"

namespace fem::solver
{
class MfemPardisoSolver final : public ILinearSolver
{
  public:
    MfemPardisoSolver() = default;

    void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) override;
};
}  // namespace fem::solver
