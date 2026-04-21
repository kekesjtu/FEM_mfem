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
#if defined(MFEM_USE_MPI)
    auto *hypre_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (!hypre_matrix)
    {
        throw std::runtime_error(
            "MfemPcgSolver requires mfem::HypreParMatrix operator (MPI parallel only).");
    }

    auto amg = std::make_unique<mfem::HypreBoomerAMG>();
    amg->SetPrintLevel(0);

    int par_print = (mfem::Mpi::WorldRank() == 0) ? print_level_ : 0;
    auto cg = std::make_unique<mfem::CGSolver>(hypre_matrix->GetComm());
    cg->SetRelTol(rel_tol_);
    cg->SetAbsTol(abs_tol_);
    cg->SetMaxIter(max_iter_);
    cg->SetPrintLevel(par_print);
    cg->SetPreconditioner(*amg);
    cg->SetOperator(*hypre_matrix);
    cg->Mult(b, x);

    if (!cg->GetConverged())
    {
        throw std::runtime_error("Parallel PCG failed to converge.");
    }
    preconditioner_ = std::move(amg);
    solver_ = std::move(cg);
#else
    (void)A;
    (void)b;
    (void)x;
    throw std::runtime_error("MfemPcgSolver requires MFEM built with MPI support.");
#endif
}

void MfemPcgSolver::Mult(const mfem::Vector &b, mfem::Vector &x)
{
    if (!solver_)
        throw std::runtime_error("MfemPcgSolver::Mult requires a prior Solve() call.");
    solver_->Mult(b, x);
    if (!solver_->GetConverged())
    {
        throw std::runtime_error("PCG failed to converge in Mult().");
    }
}
}  // namespace fem::solver
