#include "fem/fe/ParallelScalarFeContext.hpp"

#include <stdexcept>

namespace fem::fe
{
ParallelScalarFeContext::ParallelScalarFeContext(const std::string &mesh_path, int order,
                                                 int uniform_refinement_levels)
{
#if defined(MFEM_USE_MPI)
    if (!mfem::Mpi::IsInitialized())
    {
        throw std::runtime_error("MPI must be initialized before creating ParallelScalarFeContext.");
    }

    serial_mesh_ = std::make_unique<mfem::Mesh>(mesh_path.c_str(), 1, 1);
    for (int i = 0; i < uniform_refinement_levels; ++i)
    {
        serial_mesh_->UniformRefinement();
    }

    par_mesh_ = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, *serial_mesh_);
    fec_ = std::make_unique<mfem::H1_FECollection>(order, par_mesh_->Dimension());
    fes_ = std::make_unique<mfem::ParFiniteElementSpace>(par_mesh_.get(), fec_.get());
#else
    (void)mesh_path;
    (void)order;
    (void)uniform_refinement_levels;
    throw std::runtime_error("ParallelScalarFeContext requires MFEM_USE_MPI.");
#endif
}

mfem::ParMesh &ParallelScalarFeContext::Mesh()
{
    return *par_mesh_;
}

mfem::ParFiniteElementSpace &ParallelScalarFeContext::Space()
{
    return *fes_;
}
}  // namespace fem::fe
