#pragma once

#include "mfem.hpp"

#include <memory>
#include <string>

namespace fem::fe
{
class ParallelScalarFeContext
{
  public:
    ParallelScalarFeContext(const std::string &mesh_path, int order, int uniform_refinement_levels);

    mfem::ParMesh &Mesh();
    mfem::ParFiniteElementSpace &Space();

  private:
    std::unique_ptr<mfem::Mesh> serial_mesh_;
    std::unique_ptr<mfem::ParMesh> par_mesh_;
    std::unique_ptr<mfem::FiniteElementCollection> fec_;
    std::unique_ptr<mfem::ParFiniteElementSpace> fes_;
};
}  // namespace fem::fe
