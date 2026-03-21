#include "fem/solver/MfemUmfpackSolver.hpp"

#include <stdexcept>

namespace fem::solver
{
void MfemUmfpackSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    auto *sparse_matrix = dynamic_cast<mfem::SparseMatrix *>(&A);
    if (!sparse_matrix)
    {
        throw std::runtime_error("MfemUmfpackSolver expects mfem::SparseMatrix operator.");
    }

#if defined(MFEM_USE_SUITESPARSE)
    mfem::UMFPackSolver umf_solver;
    umf_solver.SetOperator(*sparse_matrix);
    umf_solver.Mult(b, x);
#else
    (void)b;
    (void)x;
    throw std::runtime_error(
        "UMFPACK solver requested, but MFEM is built without SuiteSparse (MFEM_USE_SUITESPARSE).");
#endif
}
}  // namespace fem::solver
