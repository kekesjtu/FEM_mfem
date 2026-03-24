#include "fem/app/FieldRunnerCommon.hpp"

#include "fem/app/FieldRunners.hpp"
#include "fem/fe/ParallelScalarFeContext.hpp"
#include "fem/fe/ScalarFeContext.hpp"
#include "fem/io/ParaviewExporter.hpp"
#include "fem/post/SolutionTextExporter.hpp"

#include "spdlog/spdlog.h"

#include <memory>
#include <stdexcept>

namespace fem::app
{
namespace
{
std::string ParallelAwareOutputPath(const std::string &path)
{
#if defined(MFEM_USE_MPI)
    if (mfem::Mpi::IsInitialized() && mfem::Mpi::WorldSize() > 1)
    {
        const int rank = mfem::Mpi::WorldRank();
        const std::string suffix = ".rank" + std::to_string(rank);
        const std::size_t dot = path.find_last_of('.');
        if (dot == std::string::npos)
        {
            return path + suffix;
        }
        return path.substr(0, dot) + suffix + path.substr(dot);
    }
#endif
    return path;
}

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
    const bool use_parallel_assembly = (config.simulation.assembly_mode == "parallel");

    std::unique_ptr<fem::fe::ScalarFeContext> serial_context;
    std::unique_ptr<fem::fe::ParallelScalarFeContext> parallel_context;
    mfem::Mesh *mesh = nullptr;

    if (use_parallel_assembly)
    {
        parallel_context = std::make_unique<fem::fe::ParallelScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh = &parallel_context->Mesh();
    }
    else
    {
        serial_context = std::make_unique<fem::fe::ScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh = &serial_context->Mesh();
    }

    detail::RequireDomainMaterialsCoverAllDomains(*mesh, config,
                                                  "thermo-mechanical thermal field");

    std::unordered_map<int, std::string> diffusion_by_domain;
    if (!config.domain_materials.empty())
    {
        diffusion_by_domain = detail::BuildDiffusionByDomainFromMaterials(config, materials);
    }

    const std::string diffusion_default = detail::ResolveDiffusionFromMaterial(config, materials);

    auto diffusion = detail::BuildPiecewiseDomainCoefficient(
        *mesh, diffusion_default, diffusion_by_domain, coupled_variable_names);
    auto source = detail::BuildPiecewiseDomainCoefficient(
        *mesh, thermal_field.source_default, thermal_field.source_by_domain,
        coupled_variable_names);

    std::unique_ptr<mfem::GridFunction> serial_temperature;
    std::unique_ptr<mfem::ParGridFunction> parallel_temperature;
    if (use_parallel_assembly)
    {
        auto *par_space = &parallel_context->Space();
        parallel_temperature = std::make_unique<mfem::ParGridFunction>(par_space);
        *parallel_temperature = 0.0;
        detail::SolveScalarFieldOnParallelContext(*parallel_context, thermal_field, *diffusion,
                                                  *source, 500, config.simulation.solver,
                                                  config.simulation.assembly_mode,
                                                  *parallel_temperature);
    }
    else
    {
        auto *serial_space = &serial_context->Space();
        serial_temperature = std::make_unique<mfem::GridFunction>(serial_space);
        *serial_temperature = 0.0;
        detail::SolveScalarFieldOnContext(*serial_context, thermal_field, *diffusion, *source,
                                          500, config.simulation.solver,
                                          config.simulation.assembly_mode, *serial_temperature);
    }
    mfem::GridFunction &temperature =
        use_parallel_assembly ? static_cast<mfem::GridFunction &>(*parallel_temperature)
                              : static_cast<mfem::GridFunction &>(*serial_temperature);

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        ParallelAwareOutputPath(config.simulation.output_dir + "/" + thermal_field.name +
                                "_solution.txt"),
        thermal_field.name, *mesh, temperature, "solution");

    fem::io::ParaviewExporter thermal_exporter(config.simulation.output_dir);
    thermal_exporter.Save(thermal_field.name, "temperature", *mesh, temperature, 0, 0.0);

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
    const bool use_parallel_assembly = (config.simulation.assembly_mode == "parallel");

    std::unique_ptr<fem::fe::ScalarFeContext> serial_context;
    std::unique_ptr<fem::fe::ParallelScalarFeContext> parallel_context;
    mfem::Mesh *mesh_ptr = nullptr;
    if (use_parallel_assembly)
    {
        parallel_context = std::make_unique<fem::fe::ParallelScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh_ptr = &parallel_context->Mesh();
    }
    else
    {
        serial_context = std::make_unique<fem::fe::ScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh_ptr = &serial_context->Mesh();
    }
    mfem::Mesh &mesh = *mesh_ptr;

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
    auto electro_source = detail::BuildPiecewiseDomainCoefficient(
        mesh, electro_field.source_default, electro_field.source_by_domain, coupled_variable_names);

    std::unique_ptr<mfem::GridFunction> serial_potential;
    std::unique_ptr<mfem::ParGridFunction> parallel_potential;
    if (use_parallel_assembly)
    {
        auto *par_space = &parallel_context->Space();
        parallel_potential = std::make_unique<mfem::ParGridFunction>(par_space);
        *parallel_potential = 0.0;
        detail::SolveScalarFieldOnParallelContext(*parallel_context, electro_field, *sigma_coeff,
                                                  *electro_source, 1000,
                                                  config.simulation.solver,
                                                  config.simulation.assembly_mode,
                                                  *parallel_potential);
    }
    else
    {
        auto *serial_space = &serial_context->Space();
        serial_potential = std::make_unique<mfem::GridFunction>(serial_space);
        *serial_potential = 0.0;
        detail::SolveScalarFieldOnContext(*serial_context, electro_field, *sigma_coeff,
                                          *electro_source, 1000, config.simulation.solver,
                                          config.simulation.assembly_mode, *serial_potential);
    }
    mfem::GridFunction &potential =
        use_parallel_assembly ? static_cast<mfem::GridFunction &>(*parallel_potential)
                              : static_cast<mfem::GridFunction &>(*serial_potential);

    // Only rank 0 prints logs to avoid duplication in MPI runs
#if defined(MFEM_USE_MPI)
    int rank = mfem::Mpi::WorldRank();
#else
    int rank = 0;
#endif

    detail::FieldStatistics electro_stats;
    if (use_parallel_assembly)
    {
        electro_stats = detail::ComputeGlobalFieldStatistics(*parallel_potential);
    }
    else
    {
        electro_stats = detail::ComputeFieldStatistics(*serial_potential);
    }

    if (rank == 0)
    {
        spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}",
                     electro_field.name, electro_stats.min_value, electro_stats.max_value,
                     electro_stats.mean_value, electro_stats.nan_count);
    }

    auto thermal_source_base = detail::BuildPiecewiseDomainCoefficient(
        mesh, thermal_field.source_default, thermal_field.source_by_domain, coupled_variable_names);
    auto joule_source_coeff =
        std::make_unique<JouleHeatCoefficient>(*sigma_coeff, potential, mesh.Dimension());
    auto thermal_source_raw =
        std::make_unique<SumScalarCoefficient>(*thermal_source_base, *joule_source_coeff);

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

    std::unique_ptr<mfem::GridFunction> serial_temperature;
    std::unique_ptr<mfem::ParGridFunction> parallel_temperature;
    if (use_parallel_assembly)
    {
        auto *par_space = &parallel_context->Space();
        parallel_temperature = std::make_unique<mfem::ParGridFunction>(par_space);
        *parallel_temperature = 0.0;
        detail::SolveScalarFieldOnParallelContext(*parallel_context, thermal_field,
                                                  *thermal_diffusion_coeff, *thermal_source_raw,
                                                  1000, config.simulation.solver,
                                                  config.simulation.assembly_mode,
                                                  *parallel_temperature);
    }
    else
    {
        auto *serial_space = &serial_context->Space();
        serial_temperature = std::make_unique<mfem::GridFunction>(serial_space);
        *serial_temperature = 0.0;
        detail::SolveScalarFieldOnContext(*serial_context, thermal_field,
                                          *thermal_diffusion_coeff, *thermal_source_raw, 1000,
                                          config.simulation.solver,
                                          config.simulation.assembly_mode, *serial_temperature);
    }
    mfem::GridFunction &temperature =
        use_parallel_assembly ? static_cast<mfem::GridFunction &>(*parallel_temperature)
                              : static_cast<mfem::GridFunction &>(*serial_temperature);

    detail::FieldStatistics thermal_stats;
    if (use_parallel_assembly)
    {
        thermal_stats = detail::ComputeGlobalFieldStatistics(*parallel_temperature);
    }
    else
    {
        thermal_stats = detail::ComputeFieldStatistics(*serial_temperature);
    }

    if (rank == 0)
    {
        spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}",
                     thermal_field.name, thermal_stats.min_value, thermal_stats.max_value,
                     thermal_stats.mean_value, thermal_stats.nan_count);
    }

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        ParallelAwareOutputPath(config.simulation.output_dir + "/" + electro_field.name +
                                "_solution.txt"),
        electro_field.name, mesh, potential, "solution");
    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        ParallelAwareOutputPath(config.simulation.output_dir + "/" + thermal_field.name +
                                "_solution.txt"),
        thermal_field.name, mesh, temperature, "solution");

    fem::io::ParaviewExporter exporter(config.simulation.output_dir);
    exporter.Save(electro_field.name, "solution", mesh, potential, 0, 0.0);
    exporter.Save(thermal_field.name, "solution", mesh, temperature, 0, 0.0);

    return 0;
}
}  // namespace fem::app
