#include "fem/solver/MfemPardisoSolver.hpp"

#include <stdexcept>

namespace fem::solver
{
void MfemPardisoSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    auto *hypre_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (!hypre_matrix)
    {
        throw std::runtime_error("MfemPardisoSolver requires mfem::HypreParMatrix operator "
                                 "(parallel CPARDISO path).");
    }

    auto cpardiso = std::make_unique<mfem::CPardisoSolver>(hypre_matrix->GetComm());
    cpardiso->SetPrintLevel(0);
    cpardiso->SetMatrixType(mfem::CPardisoSolver::REAL_STRUCTURE_SYMMETRIC);
    cpardiso->SetOperator(*hypre_matrix);
    cpardiso->Mult(b, x);
    solver_ = std::move(cpardiso);
}

void MfemPardisoSolver::Mult(const mfem::Vector &b, mfem::Vector &x)
{
    if (!solver_)
        throw std::runtime_error("MfemPardisoSolver::Mult requires a prior Solve() call.");
    solver_->Mult(b, x);
}
}  // namespace fem::solver
