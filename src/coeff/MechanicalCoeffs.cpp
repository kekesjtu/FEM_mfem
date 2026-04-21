#include "fem/coeff/MechanicalCoeffs.hpp"

#include "fem/coeff/CoefficientManager.hpp"
#include "fem/frontend/Config.hpp"

#include <stdexcept>

namespace fem::coeff
{

void MechanicalCoeffs::BuildCoeffs(frontend::ProjectConfig &config)
{
    if (lambda_ && mu_ && body_force_ && alpha_ && reference_temperature_coeff_)
        return;

    auto &mesh = config.fe.GetMesh();
    reference_temperature_value_ = config.mechanical_field.reference_temperature;

    lambda_ = std::make_unique<coeff::PiecewiseConstantCoefficient>("lambda_lame", config.materials,
                                                                    mesh);
    mu_ = std::make_unique<coeff::PiecewiseConstantCoefficient>("mu_lame", config.materials, mesh);

    const auto &field_config = config.mechanical_field;
    const int dim = mesh.Dimension();

    if (!field_config.domain_to_body_force.empty())
    {
        int max_attr = mesh.attributes.Max();
        auto pw_vec = std::make_unique<mfem::VectorArrayCoefficient>(dim);
        for (int d = 0; d < dim; ++d)
        {
            mfem::Vector comp_vals(max_attr);
            double default_val = (d < static_cast<int>(field_config.body_force_default.size()))
                                     ? field_config.body_force_default[d]
                                     : 0.0;
            for (int a = 0; a < max_attr; ++a)
            {
                auto it = field_config.domain_to_body_force.find(a + 1);
                if (it != field_config.domain_to_body_force.end() &&
                    d < static_cast<int>(it->second.size()))
                    comp_vals(a) = it->second[d];
                else
                    comp_vals(a) = default_val;
            }
            pw_vec->Set(d, new mfem::PWConstCoefficient(comp_vals));
        }
        body_force_ = std::move(pw_vec);
    }
    else
    {
        mfem::Vector body_force_vec(dim);
        body_force_vec = 0.0;
        for (int d = 0; d < dim && d < static_cast<int>(field_config.body_force_default.size());
             ++d)
            body_force_vec(d) = field_config.body_force_default[d];
        body_force_ = std::make_unique<mfem::VectorConstantCoefficient>(body_force_vec);
    }

    if (!alpha_)
    {
        alpha_ = std::make_unique<coeff::PiecewiseConstantCoefficient>(
            "thermal_expansion_secant", config.materials, config.fe.GetMesh());
    }

    if (!reference_temperature_coeff_)
    {
        reference_temperature_coeff_ =
            std::make_unique<mfem::ConstantCoefficient>(reference_temperature_value_);
    }

    // default to reference temperature if no field provided; SetTemperatureField can be called
    // later to update
    if (!temperature_)
    {
        temperature_ = reference_temperature_coeff_.get();
    }
}

void MechanicalCoeffs::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    if (!temperature_gf)
        throw std::invalid_argument("MechanicalCoeffs::SetTemperatureField does not accept null");

    if (!field_temperature_coeff_)
    {
        field_temperature_coeff_ = std::make_unique<mfem::GridFunctionCoefficient>(temperature_gf);
    }
    else
    {
        field_temperature_coeff_->SetGridFunction(temperature_gf);
    }

    temperature_ = field_temperature_coeff_.get();
}

}  // namespace fem::coeff
