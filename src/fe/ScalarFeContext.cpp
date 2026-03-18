#include "fem/fe/ScalarFeContext.hpp"

namespace fem::fe
{
ScalarFeContext::ScalarFeContext(const std::string &mesh_path, int order,
                                 int uniform_refinement_levels)
{
    mesh_ = std::make_unique<mfem::Mesh>(mesh_path.c_str(), 1, 1);
    for (int i = 0; i < uniform_refinement_levels; ++i)
    {
        mesh_->UniformRefinement();
    }

    fec_ = std::make_unique<mfem::H1_FECollection>(order, mesh_->Dimension());
    fes_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());
}

mfem::Mesh &ScalarFeContext::Mesh()
{
    return *mesh_;
}

mfem::FiniteElementSpace &ScalarFeContext::Space()
{
    return *fes_;
}
}  // namespace fem::fe
