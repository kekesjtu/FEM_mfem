#include "fem/solver/MfemUmfpackSolver.hpp"

#include <stdexcept>

namespace fem::solver
{

void MfemUmfpackSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    (void)A;
    (void)b;
    (void)x;
    throw std::runtime_error(
        "MfemUmfpackSolver: UMFPACK is serial-only and not supported in parallel configuration.");
}

void MfemUmfpackSolver::Mult(const mfem::Vector &b, mfem::Vector &x)
{
    (void)b;
    (void)x;
    throw std::runtime_error(
        "MfemUmfpackSolver: UMFPACK is serial-only and not supported in parallel configuration.");
}
}  // namespace fem::solver
