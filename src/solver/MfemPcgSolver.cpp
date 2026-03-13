#include "fem/solver/MfemPcgSolver.hpp"

#include <stdexcept>

namespace fem::solver
{
MfemPcgSolver::MfemPcgSolver(double rel_tol, double abs_tol, int max_iter, int print_level)
    : rel_tol_(rel_tol), abs_tol_(abs_tol), max_iter_(max_iter), print_level_(print_level)
{
}

void MfemPcgSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    auto *sparse_matrix = dynamic_cast<mfem::SparseMatrix *>(&A);
    if (!sparse_matrix)
    {
        throw std::runtime_error("MfemPcgSolver expects mfem::SparseMatrix operator.");
    }

    mfem::GSSmoother preconditioner(*sparse_matrix);
    mfem::PCG(*sparse_matrix, preconditioner, b, x, print_level_, max_iter_, rel_tol_, abs_tol_);
}
}  // namespace fem::solver
