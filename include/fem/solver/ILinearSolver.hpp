#pragma once

#include "mfem.hpp"

namespace fem::solver
{
class ILinearSolver
{
  public:
    virtual ~ILinearSolver() = default;
    virtual void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) = 0;
};
}  // namespace fem::solver
