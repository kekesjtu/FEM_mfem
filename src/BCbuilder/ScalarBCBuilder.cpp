#include "fem/BCbuilder/ScalarBCBuilder.hpp"

#include "fem/frontend/Config.hpp"

namespace fem::BCbuilder
{

ScalarBCData ScalarBCBuilder::BuildBC(const frontend::ScalarFieldConfig &field_config,
                                      const mfem::ParMesh &mesh,
                                      mfem::ParFiniteElementSpace &fespace)
{
    ScalarBCData bc;

    const int num_bdr = mesh.bdr_attributes.Max();

    bc.essential_bdr.SetSize(num_bdr);
    bc.essential_bdr = 0;
    bc.dirichlet_coeffs.reserve(field_config.dirichlet_bcs.size());
    bc.dirichlet_markers.reserve(field_config.dirichlet_bcs.size());
    for (const auto &dbc : field_config.dirichlet_bcs)
    {
        auto marker = BuildBoundaryMarker(num_bdr, dbc.bdr_attributes, &bc.essential_bdr);
        bc.dirichlet_coeffs.emplace_back(dbc.value);
        bc.dirichlet_markers.push_back(std::move(marker));
    }

    bc.robin_l_coeffs.reserve(field_config.robin_bcs.size());
    bc.robin_q_coeffs.reserve(field_config.robin_bcs.size());
    bc.robin_markers.reserve(field_config.robin_bcs.size());
    for (const auto &rbc : field_config.robin_bcs)
    {
        auto marker = BuildBoundaryMarker(num_bdr, rbc.bdr_attributes);
        bc.robin_l_coeffs.emplace_back(rbc.l);
        bc.robin_q_coeffs.emplace_back(rbc.q);
        bc.robin_markers.push_back(std::move(marker));
    }

    fespace.GetEssentialTrueDofs(bc.essential_bdr, bc.essential_tdofs);
    return bc;
}

}  // namespace fem::BCbuilder
