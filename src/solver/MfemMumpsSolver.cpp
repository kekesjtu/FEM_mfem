#include "fem/solver/MfemMumpsSolver.hpp"

#include <stdexcept>

#include <memory>

namespace fem::solver
{
// mumps is only
void MfemMumpsSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    auto *hypre_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (!hypre_matrix)
    {
        throw std::runtime_error(
            "MfemMumpsSolver requires mfem::HypreParMatrix operator (MPI parallel only).");
    }

    auto mumps = std::make_unique<mfem::MUMPSSolver>(hypre_matrix->GetComm());
    mumps->SetPrintLevel(0);
    mumps->SetMatrixSymType(mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE);
    mumps->SetOperator(*hypre_matrix);
    mumps->Mult(b, x);
    solver_ = std::move(mumps);
}

void MfemMumpsSolver::Mult(const mfem::Vector &b, mfem::Vector &x)
{
    if (!solver_)
        throw std::runtime_error("MfemMumpsSolver::Mult requires a prior Solve() call.");
    solver_->Mult(b, x);
}
}  // namespace fem::solver
