#include "fem/coeff/ElectrostaticCoeffs.hpp"

#include "fem/coeff/CoefficientManager.hpp"
#include "fem/frontend/Config.hpp"

namespace fem::coeff
{

void ElectrostaticCoeffs::BuildCoeffs(frontend::ProjectConfig &config,
                                      mfem::GridFunction *const *temperature_gf,
                                      mfem::Coefficient *const *sigma)
{
    temperature_gf_ = temperature_gf;
    sigma_ = sigma;

    bc_expressions_.clear();
    bc_expressions_.reserve(config.electric_field.dirichlet_bcs.size());
    for (const auto &dbc : config.electric_field.dirichlet_bcs)
    {
        if (!dbc.value_expr.empty())
            bc_expressions_.emplace_back(dbc.value_expr);
        else
            bc_expressions_.emplace_back(std::to_string(dbc.value));
    }

    if (!owned_sigma_)
    {
        auto local_sigma = std::make_unique<coeff::ExpressionCoefficient>(
            "electrical_conductivity", config.materials, config.fe.GetMesh(), temperature_gf_);
        local_sigma->SetReferenceTemperature(config.electric_field.reference_temperature);
        owned_sigma_ = std::move(local_sigma);
    }
}

void ElectrostaticCoeffs::UpdateBoundaryConditions(const frontend::ScalarFieldConfig &field_config,
                                                   double t, BCbuilder::ScalarBCData &bc) const
{
    for (size_t i = 0; i < field_config.dirichlet_bcs.size() && i < bc_expressions_.size() &&
                       i < bc.dirichlet_coeffs.size();
         ++i)
    {
        if (!field_config.dirichlet_bcs[i].value_expr.empty())
        {
            double val = bc_expressions_[i].Evaluate({0.0, 0.0, 0.0, t, 0.0});
            bc.dirichlet_coeffs[i].constant = val;
        }
    }
}

}  // namespace fem::coeff
