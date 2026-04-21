#include "fem/coeff/ThermalCoeffs.hpp"

#include "fem/frontend/Config.hpp"

#include <string>

namespace fem::coeff
{

void ThermalCoeffs::BuildCoeffs(frontend::ProjectConfig &config,
                                mfem::GridFunction *const *voltage_gf,
                                mfem::Coefficient *const *sigma)
{
    voltage_gf_ = voltage_gf;
    sigma_ = sigma;

    if (!k_ || !rho_cp_)
    {
        auto &mesh = config.fe.GetMesh();
        k_ = std::make_unique<coeff::PiecewiseConstantCoefficient>("diffusion", config.materials,
                                                                   mesh);
        rho_cp_ =
            std::make_unique<coeff::PiecewiseConstantCoefficient>("rho_cp", config.materials, mesh);
    }

    if (!base_source_)
    {
        const auto &field_config = config.thermal_field;
        if (!field_config.domain_to_source.empty())
        {
            int max_attr = config.fe.GetMesh().attributes.Max();
            auto pw = std::make_unique<mfem::PWConstCoefficient>(max_attr);
            double default_source = std::stod(field_config.source_default);
            for (int a = 1; a <= max_attr; ++a)
            {
                auto it = field_config.domain_to_source.find(a);
                pw->operator()(a) = (it != field_config.domain_to_source.end())
                                        ? std::stod(it->second)
                                        : default_source;
            }
            base_source_ = std::move(pw);
        }
        else
        {
            base_source_ =
                std::make_unique<mfem::ConstantCoefficient>(std::stod(field_config.source_default));
        }
    }

    if (!joule_source_ && config.HasElectricField())
        joule_source_ = std::make_unique<coeff::JouleHeatingCoefficient>(sigma_, voltage_gf_);

    if (!total_source_)
    {
        total_source_ = std::make_unique<coeff::CombinedHeatSourceCoefficient>(base_source_.get(),
                                                                               joule_source_.get());
    }
}

}  // namespace fem::coeff
