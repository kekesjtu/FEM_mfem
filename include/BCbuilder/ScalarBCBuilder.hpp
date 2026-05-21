#pragma once

#include "mfem.hpp"

#include <vector>

namespace fem::frontend
{
struct ScalarFieldConfig;
}

namespace fem::BCbuilder
{

struct ScalarBCData
{
    mfem::Array<int> essential_bdr;
    mfem::Array<int> essential_tdofs;
    std::vector<mfem::ConstantCoefficient> dirichlet_coeffs;
    std::vector<mfem::Array<int>> dirichlet_markers;
    std::vector<mfem::ConstantCoefficient> robin_l_coeffs;
    std::vector<mfem::ConstantCoefficient> robin_q_coeffs;
    std::vector<mfem::Array<int>> robin_markers;
};

inline mfem::Array<int> BuildBoundaryMarker(int num_bdr_attributes,
                                            const std::vector<int> &bdr_attributes,
                                            mfem::Array<int> *essential_bdr = nullptr)
{
    mfem::Array<int> marker(num_bdr_attributes);
    marker = 0;
    for (int bdr_attr : bdr_attributes)
    {
        if (bdr_attr >= 1 && bdr_attr <= num_bdr_attributes)
        {
            marker[bdr_attr - 1] = 1;
            if (essential_bdr)
                (*essential_bdr)[bdr_attr - 1] = 1;
        }
    }
    return marker;
}

class ScalarBCBuilder
{
  public:
    static ScalarBCData BuildBC(const frontend::ScalarFieldConfig &field_config,
                                const mfem::ParMesh &mesh, mfem::ParFiniteElementSpace &fespace);
};

}  // namespace fem::BCbuilder
