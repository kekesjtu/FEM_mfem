#include "fem/solver/MfemAmgSolver.hpp"

#include <stdexcept>

namespace fem::solver
{
MfemAmgSolver::MfemAmgSolver(double rel_tol, double abs_tol, int max_iter, int print_level)
    : rel_tol_(rel_tol), abs_tol_(abs_tol), max_iter_(max_iter), print_level_(print_level)
{
}

void MfemAmgSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
#if defined(MFEM_USE_MPI)
    auto *hypre_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (!hypre_matrix)
    {
        throw std::runtime_error(
            "MfemAmgSolver expects mfem::HypreParMatrix operator. "
            "Current assembly provides non-hypre matrix; use solver='pcg' or switch to parallel hypre assembly.");
    }

    mfem::HypreBoomerAMG amg;
    amg.SetPrintLevel(print_level_ > 0 ? 1 : 0);
    amg.SetOperator(*hypre_matrix);

    mfem::HyprePCG solver(hypre_matrix->GetComm());
    solver.SetTol(rel_tol_);
    solver.SetAbsTol(abs_tol_);
    solver.SetMaxIter(max_iter_);
    solver.SetPrintLevel(print_level_);
    solver.SetPreconditioner(amg);
    solver.SetOperator(*hypre_matrix);
    solver.Mult(b, x);
#else
    (void)A;
    (void)b;
    (void)x;
    throw std::runtime_error(
        "AMG solver requested, but current MFEM build has no MPI/Hypre support. "
        "Rebuild MFEM with MPI+Hypre, or use solver='pcg' in current serial setup.");
#endif
}
}  // namespace fem::solver
