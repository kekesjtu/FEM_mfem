#pragma once

#include "mfem.hpp"

#include <vector>

namespace fem::frontend
{
struct ProjectConfig;
}

namespace fem::BCbuilder
{

struct MechanicalBCData
{
    mfem::Array<int> essential_bdr;
    mfem::Array<int> essential_tdofs;
    std::vector<mfem::VectorConstantCoefficient> disp_coeffs;
    std::vector<mfem::Array<int>> disp_markers;
    std::vector<mfem::ConstantCoefficient> penalty_coeffs;
    std::vector<mfem::Array<int>> nd_markers;
    std::vector<mfem::ConstantCoefficient> pressure_coeffs;
    std::vector<mfem::Array<int>> pressure_markers;
};

class MechanicalBCBuilder
{
  public:
    static MechanicalBCData BuildBC(const frontend::ProjectConfig &config);
};

}  // namespace fem::BCbuilder