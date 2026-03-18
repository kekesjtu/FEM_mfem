#include "fem/app/FieldRunnerCommon.hpp"

#include "fem/app/FieldRunners.hpp"
#include "fem/fe/ScalarFeContext.hpp"
#include "fem/io/ParaviewExporter.hpp"
#include "fem/post/SolutionTextExporter.hpp"

#include "spdlog/spdlog.h"

#include <memory>

namespace fem::app
{
namespace
{
class JouleHeatCoefficient : public mfem::Coefficient
{
  public:
    JouleHeatCoefficient(mfem::Coefficient &sigma, const mfem::GridFunction &potential, int dim)
        : sigma_(sigma), potential_(potential), dim_(dim)
    {
    }

    double Eval(mfem::ElementTransformation &trans, const mfem::IntegrationPoint &ip) override
    {
        trans.SetIntPoint(&ip);
        mfem::Vector grad(dim_);
        potential_.GetGradient(trans, grad);
        const double sigma = sigma_.Eval(trans, ip);
        return sigma * (grad * grad);
    }

  private:
    mfem::Coefficient &sigma_;
    const mfem::GridFunction &potential_;
    int dim_;
};

class SumScalarCoefficient : public mfem::Coefficient
{
  public:
    SumScalarCoefficient(mfem::Coefficient &lhs, mfem::Coefficient &rhs) : lhs_(lhs), rhs_(rhs)
    {
    }

    double Eval(mfem::ElementTransformation &trans, const mfem::IntegrationPoint &ip) override
    {
        return lhs_.Eval(trans, ip) + rhs_.Eval(trans, ip);
    }

  private:
    mfem::Coefficient &lhs_;
    mfem::Coefficient &rhs_;
};
}  // namespace

int RunThermoMechanicalCoupled(const fem::frontend::ProjectConfig &config,
                               const fem::frontend::FieldConfig &thermal_field,
                               const fem::frontend::FieldConfig &mechanical_field,
                               const std::vector<std::string> &coupled_variable_names,
                               fem::material::MaterialDatabase &materials)
{
    fem::fe::ScalarFeContext thermal_context(config.simulation.mesh_path, config.simulation.order,
                                             config.simulation.uniform_refinement_levels);

    detail::RequireDomainMaterialsCoverAllDomains(thermal_context.Mesh(), config,
                                                  "thermo-mechanical thermal field");

    std::unordered_map<int, std::string> diffusion_by_domain;
    if (!config.domain_materials.empty())
    {
        diffusion_by_domain = detail::BuildDiffusionByDomainFromMaterials(config, materials);
    }

    const std::string diffusion_default = detail::ResolveDiffusionFromMaterial(config, materials);

    auto diffusion = detail::BuildPiecewiseDomainCoefficient(
        thermal_context.Mesh(), diffusion_default, diffusion_by_domain, coupled_variable_names);
    auto source = detail::BuildPiecewiseDomainCoefficient(
        thermal_context.Mesh(), thermal_field.source_default, thermal_field.source_by_domain,
        coupled_variable_names);

    mfem::GridFunction temperature(&thermal_context.Space());
    temperature = 0.0;

    detail::SolveScalarFieldOnContext(thermal_context, thermal_field, *diffusion, *source, 1000,
                                      temperature);

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        config.simulation.output_dir + "/" + thermal_field.name + "_solution.txt",
        thermal_field.name, thermal_context.Mesh(), temperature, "solution");

    fem::io::ParaviewExporter thermal_exporter(config.simulation.output_dir);
    thermal_exporter.Save(thermal_field.name, "temperature", thermal_context.Mesh(), temperature, 0,
                          0.0);

    mfem::GridFunctionCoefficient temperature_field_coefficient(&temperature);

    return detail::RunMechanicalFieldInternal(config, mechanical_field, coupled_variable_names,
                                              materials, &temperature_field_coefficient);
}

int RunElectroThermalCoupled(const fem::frontend::ProjectConfig &config,
                             const fem::frontend::FieldConfig &electro_field,
                             const fem::frontend::FieldConfig &thermal_field,
                             const std::vector<std::string> &coupled_variable_names,
                             fem::material::MaterialDatabase &materials)
{
    fem::fe::ScalarFeContext electro_context(config.simulation.mesh_path, config.simulation.order,
                                             config.simulation.uniform_refinement_levels);
    mfem::Mesh &mesh = electro_context.Mesh();

    detail::RequireDomainMaterialsCoverAllDomains(mesh, config, "electro-thermal coupled fields");

    std::unordered_map<int, std::string> sigma_by_domain_expr;
    if (!config.domain_materials.empty())
    {
        sigma_by_domain_expr = detail::BuildPropertyByDomainFromMaterials(
            config, materials, "electrical_conductivity");
        if (sigma_by_domain_expr.empty())
        {
            sigma_by_domain_expr =
                detail::BuildPropertyByDomainFromMaterials(config, materials, "diffusion");
            if (!sigma_by_domain_expr.empty())
            {
                spdlog::warn("No 'electrical_conductivity' found in domain materials; fallback to "
                             "'diffusion' for electrostatic conductivity.");
            }
        }
    }

    const std::string sigma_default =
        detail::ResolveElectricalConductivityFromMaterial(config, materials);

    auto sigma_coeff = detail::BuildPiecewiseDomainCoefficient(
        mesh, sigma_default, sigma_by_domain_expr, coupled_variable_names);
    const mfem::Vector sigma_values = detail::BuildPiecewiseDomainValues(
        mesh, sigma_default, sigma_by_domain_expr, coupled_variable_names);

    auto electro_source = detail::BuildPiecewiseDomainCoefficient(
        mesh, electro_field.source_default, electro_field.source_by_domain, coupled_variable_names);

    mfem::GridFunction potential(&electro_context.Space());
    potential = 0.0;

    detail::SolveScalarFieldOnContext(electro_context, electro_field, *sigma_coeff, *electro_source,
                                      1000, potential);

    const auto electro_stats = detail::ComputeFieldStatistics(potential);
    spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}",
                 electro_field.name, electro_stats.min_value, electro_stats.max_value,
                 electro_stats.mean_value, electro_stats.nan_count);

    auto thermal_source_base = detail::BuildPiecewiseDomainCoefficient(
        mesh, thermal_field.source_default, thermal_field.source_by_domain, coupled_variable_names);
    auto joule_source_coeff =
        std::make_unique<JouleHeatCoefficient>(*sigma_coeff, potential, mesh.Dimension());
    auto thermal_source_coeff =
        std::make_unique<SumScalarCoefficient>(*thermal_source_base, *joule_source_coeff);

    const mfem::Vector avg_grad_sq = detail::BuildDomainAverageGradSquare(mesh, potential);
    mfem::Vector joule_q(avg_grad_sq.Size());
    for (int i = 0; i < joule_q.Size(); ++i)
    {
        joule_q(i) = sigma_values(i) * avg_grad_sq(i);
    }
    detail::DebugLogDomainJouleSource(thermal_field.name, mesh, sigma_values, avg_grad_sq, joule_q);

    std::unordered_map<int, std::string> thermal_diffusion_by_domain;
    if (!config.domain_materials.empty())
    {
        thermal_diffusion_by_domain =
            detail::BuildDiffusionByDomainFromMaterials(config, materials);
    }
    const std::string thermal_diffusion_default =
        detail::ResolveDiffusionFromMaterial(config, materials);
    auto thermal_diffusion_coeff = detail::BuildPiecewiseDomainCoefficient(
        mesh, thermal_diffusion_default, thermal_diffusion_by_domain, coupled_variable_names);

    mfem::GridFunction temperature(&electro_context.Space());
    temperature = 0.0;

    detail::SolveScalarFieldOnContext(electro_context, thermal_field, *thermal_diffusion_coeff,
                                      *thermal_source_coeff, 1000, temperature);

    const auto thermal_stats = detail::ComputeFieldStatistics(temperature);
    spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}",
                 thermal_field.name, thermal_stats.min_value, thermal_stats.max_value,
                 thermal_stats.mean_value, thermal_stats.nan_count);

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        config.simulation.output_dir + "/" + electro_field.name + "_solution.txt",
        electro_field.name, mesh, potential, "solution");
    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        config.simulation.output_dir + "/" + thermal_field.name + "_solution.txt",
        thermal_field.name, mesh, temperature, "solution");

    fem::io::ParaviewExporter exporter(config.simulation.output_dir);
    exporter.Save(electro_field.name, "solution", mesh, potential, 0, 0.0);
    exporter.Save(thermal_field.name, "solution", mesh, temperature, 0, 0.0);

    return 0;
}
}  // namespace fem::app
