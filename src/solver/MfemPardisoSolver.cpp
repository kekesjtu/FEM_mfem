#include "fem/solver/MfemPardisoSolver.hpp"

#include <stdexcept>

namespace fem::solver
{

// Extract upper triangular part of a SparseMatrix for PARDISO SPD mode (mtype=2).
// PARDISO symmetric modes require only upper triangular CSR input.
static mfem::SparseMatrix* ExtractUpperTriangular(mfem::SparseMatrix &full)
{
    const int n = full.Size();
    const int *I = full.HostReadI();
    const int *J = full.HostReadJ();
    const double *V = full.HostReadData();

    int ut_nnz = 0;
    for (int i = 0; i < n; i++)
        for (int jj = I[i]; jj < I[i+1]; jj++)
            if (J[jj] >= i) ut_nnz++;

    int *ut_I = new int[n + 1];
    int *ut_J = new int[ut_nnz];
    double *ut_V = new double[ut_nnz];

    ut_I[0] = 0;
    int pos = 0;
    for (int i = 0; i < n; i++)
    {
        for (int jj = I[i]; jj < I[i+1]; jj++)
        {
            if (J[jj] >= i)
            {
                ut_J[pos] = J[jj];
                ut_V[pos] = V[jj];
                pos++;
            }
        }
        ut_I[i+1] = pos;
    }

    return new mfem::SparseMatrix(ut_I, ut_J, ut_V, n, n);
}

void MfemPardisoSolver::Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x)
{
    // Try parallel path first (HypreParMatrix → Cluster PARDISO)
#if defined(MFEM_USE_MKL_CPARDISO) && defined(MFEM_USE_MPI)
    auto *par_matrix = dynamic_cast<mfem::HypreParMatrix *>(&A);
    if (par_matrix)
    {
        mfem::CPardisoSolver cpardiso(par_matrix->GetComm());
        cpardiso.SetPrintLevel(0);
        cpardiso.SetMatrixType(mfem::CPardisoSolver::REAL_STRUCTURE_SYMMETRIC);
        cpardiso.SetOperator(*par_matrix);
        cpardiso.Mult(b, x);
        return;
    }
#endif

    // Serial path: extract upper triangular + PARDISO SPD mode
#if defined(MFEM_USE_MKL_PARDISO)
    auto *sparse_matrix = dynamic_cast<mfem::SparseMatrix *>(&A);
    if (sparse_matrix)
    {
        auto *ut_mat = ExtractUpperTriangular(*sparse_matrix);
        mfem::PardisoSolver pardiso;
        pardiso.SetPrintLevel(0);
        pardiso.SetMatrixType(mfem::PardisoSolver::REAL_SYMMETRIC_POSITIVE_DEFINITE);
        pardiso.SetOperator(*ut_mat);
        pardiso.Mult(b, x);
        delete ut_mat;
        return;
    }
#endif

    throw std::runtime_error("MfemPardisoSolver: unsupported operator type or MFEM built without "
                             "MKL PARDISO/CPARDISO support.");
}
}  // namespace fem::solver
