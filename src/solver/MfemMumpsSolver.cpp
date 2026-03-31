#include "fem/solver/MfemMumpsSolver.hpp"

#include <stdexcept>

#include <memory>

namespace fem::solver
{
void MfemMumpsSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
#if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
    auto *par_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (par_matrix)
    {
        mfem::MUMPSSolver mumps(par_matrix->GetComm());
        mumps.SetPrintLevel(0);
        mumps.SetMatrixSymType(mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE);
        mumps.SetOperator(*par_matrix);
        mumps.Mult(b, x);
        return;
    }

    auto *sparse_matrix = dynamic_cast<mfem::SparseMatrix *>(&A);
    if (sparse_matrix)
    {
        std::unique_ptr<mfem::SparseMatrix> local_matrix(new mfem::SparseMatrix(*sparse_matrix));
        HYPRE_BigInt row_starts[2] = {0, static_cast<HYPRE_BigInt>(local_matrix->Height())};
        mfem::HypreParMatrix par_self(MPI_COMM_SELF,
                                      static_cast<HYPRE_BigInt>(local_matrix->Height()), row_starts,
                                      local_matrix.get());

        mfem::MUMPSSolver mumps(MPI_COMM_SELF);
        mumps.SetPrintLevel(0);
        mumps.SetMatrixSymType(mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE);
        mumps.SetOperator(par_self);
        mumps.Mult(b, x);
        return;
    }

    throw std::runtime_error(
        "MfemMumpsSolver expects mfem::HypreParMatrix or mfem::SparseMatrix operator.");
#else
    (void)A;
    (void)b;
    (void)x;
    throw std::runtime_error(
        "MUMPS solver requested, but MFEM is built without MUMPS+MPI support.");
#endif
}
}  // namespace fem::solver
