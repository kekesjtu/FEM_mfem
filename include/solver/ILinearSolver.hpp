#pragma once

#include "mfem.hpp"

namespace fem::solver
{
class ILinearSolver
{
  public:
    virtual ~ILinearSolver() = default;
    virtual void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) = 0;

    /// Back-substitute with the operator from the last Solve() call.
    /// Direct solvers reuse the cached factorization; iterative solvers re-solve.
    virtual void Mult(const mfem::Vector &b, mfem::Vector &x) = 0;
};
}  // namespace fem::solver
