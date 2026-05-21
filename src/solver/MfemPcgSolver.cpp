#include "solver/MfemPcgSolver.hpp"

#include "log/Logger.hpp"

#include <stdexcept>

namespace fem::solver
{
MfemPcgSolver::MfemPcgSolver(double rel_tol, double abs_tol, int max_iter, int print_level,
                             mfem::ParFiniteElementSpace *elasticity_fespace)
    : rel_tol_(rel_tol),
      abs_tol_(abs_tol),
      max_iter_(max_iter),
      print_level_(print_level),
      elasticity_fespace_(elasticity_fespace)
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
    cg->iterative_mode = elasticity_fespace_ != nullptr;
    cg->SetPreconditioner(*amg);
    cg->SetOperator(*hypre_matrix);
    if (elasticity_fespace_)
    {
        const int dim = elasticity_fespace_->GetParMesh()->Dimension();
        amg->SetSystemsOptions(dim,
                               elasticity_fespace_->GetOrdering() == mfem::Ordering::byNODES);
    }
    cg->Mult(b, x);

    if (!cg->GetConverged())
    {
        throw std::runtime_error("Parallel PCG failed to converge: iterations=" +
                                 std::to_string(cg->GetNumIterations()) +
                                 ", final_rel_norm=" +
                                 std::to_string(cg->GetFinalRelNorm()));
    }
    if (print_level_ >= 2 && mfem::Mpi::WorldRank() == 0)
    {
        auto logger = fem::log::Get();
        logger->info("PCG converged: iterations={}, final_rel_norm={:.3e}",
                     cg->GetNumIterations(), cg->GetFinalRelNorm());
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
        throw std::runtime_error("PCG failed to converge in Mult(): iterations=" +
                                 std::to_string(solver_->GetNumIterations()) +
                                 ", final_rel_norm=" +
                                 std::to_string(solver_->GetFinalRelNorm()));
    }
    if (print_level_ >= 2 && mfem::Mpi::WorldRank() == 0)
    {
        auto logger = fem::log::Get();
        logger->info("PCG converged: iterations={}, final_rel_norm={:.3e}",
                     solver_->GetNumIterations(), solver_->GetFinalRelNorm());
    }
}
}  // namespace fem::solver
