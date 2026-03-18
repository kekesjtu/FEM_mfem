#pragma once

#include "mfem.hpp"

#include <memory>
#include <string>

namespace fem::fe
{
class ScalarFeContext
{
  public:
    ScalarFeContext(const std::string &mesh_path, int order, int uniform_refinement_levels);

    mfem::Mesh &Mesh();
    mfem::FiniteElementSpace &Space();

  private:
    std::unique_ptr<mfem::Mesh> mesh_;
    std::unique_ptr<mfem::FiniteElementCollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
};
}  // namespace fem::fe
