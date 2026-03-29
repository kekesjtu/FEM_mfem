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
    if (sparse_matrix)
    {
        mfem::DSmoother preconditioner(*sparse_matrix);

        mfem::CGSolver solver;
        solver.SetRelTol(rel_tol_);
        solver.SetAbsTol(abs_tol_);
        solver.SetMaxIter(max_iter_);
        solver.SetPrintLevel(print_level_);
        solver.SetPreconditioner(preconditioner);
        solver.SetOperator(*sparse_matrix);
        solver.Mult(b, x);
        if (!solver.GetConverged())
        {
            throw std::runtime_error("PCG failed to converge on serial SparseMatrix.");
        }
        return;
    }

#if defined(MFEM_USE_MPI)
    auto *hypre_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (hypre_matrix)
    {
        mfem::HypreBoomerAMG amg;
        amg.SetPrintLevel(0);

        int par_print = (mfem::Mpi::WorldRank() == 0) ? print_level_ : 0;
        mfem::CGSolver solver(hypre_matrix->GetComm());
        solver.SetRelTol(rel_tol_);
        solver.SetAbsTol(abs_tol_);
        solver.SetMaxIter(max_iter_);
        solver.SetPrintLevel(par_print);
        solver.SetPreconditioner(amg);
        solver.SetOperator(*hypre_matrix);
        solver.Mult(b, x);

        if (!solver.GetConverged())
        {
            throw std::runtime_error("Parallel PCG failed to converge.");
        }
        return;
    }
#endif

    throw std::runtime_error(
        "MfemPcgSolver expects mfem::SparseMatrix (serial) or mfem::HypreParMatrix (parallel).");
}
}  // namespace fem::solver
