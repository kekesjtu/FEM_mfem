#include "fem/BCbuilder/MechanicalBCBuilder.hpp"

#include "fem/BCbuilder/ScalarBCBuilder.hpp"

#include "fem/frontend/Config.hpp"

namespace fem::BCbuilder
{

MechanicalBCData MechanicalBCBuilder::BuildBC(const frontend::ProjectConfig &config)
{
    MechanicalBCData bc;

    auto &mesh = config.fe.GetMesh();
    const auto &field_config = config.mechanical_field;
    const int num_bdr = mesh.bdr_attributes.Max();
    const int dim = mesh.Dimension();

    bc.essential_bdr.SetSize(num_bdr);
    bc.essential_bdr = 0;
    bc.disp_coeffs.reserve(field_config.displacement_bcs.size());
    bc.disp_markers.reserve(field_config.displacement_bcs.size());
    for (const auto &dbc : field_config.displacement_bcs)
    {
        auto marker = BuildBoundaryMarker(num_bdr, dbc.bdr_attributes, &bc.essential_bdr);

        mfem::Vector val(dim);
        val = 0.0;
        for (int d = 0; d < dim && d < static_cast<int>(dbc.value.size()); ++d)
            val(d) = dbc.value[d];

        bc.disp_coeffs.emplace_back(val);
        bc.disp_markers.push_back(std::move(marker));
    }

    bc.penalty_coeffs.reserve(field_config.normal_displacement_bcs.size());
    bc.nd_markers.reserve(field_config.normal_displacement_bcs.size());
    for (const auto &nbc : field_config.normal_displacement_bcs)
    {
        auto marker = BuildBoundaryMarker(num_bdr, nbc.bdr_attributes);
        bc.penalty_coeffs.emplace_back(nbc.penalty);
        bc.nd_markers.push_back(std::move(marker));
    }

    bc.pressure_coeffs.reserve(field_config.pressure_bcs.size());
    bc.pressure_markers.reserve(field_config.pressure_bcs.size());
    for (const auto &pbc : field_config.pressure_bcs)
    {
        auto marker = BuildBoundaryMarker(num_bdr, pbc.bdr_attributes);
        bc.pressure_coeffs.emplace_back(pbc.value);
        bc.pressure_markers.push_back(std::move(marker));
    }

    config.fe.GetVectorFESpace().GetEssentialTrueDofs(bc.essential_bdr, bc.essential_tdofs);

    return bc;
}

}  // namespace fem::BCbuilder