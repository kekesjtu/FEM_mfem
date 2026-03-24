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
}  // namespace

int RunScalarField(const fem::frontend::ProjectConfig &config,
                   const fem::frontend::FieldConfig &field,
                   const std::vector<std::string> &coupled_variable_names,
                   fem::material::MaterialDatabase &materials)
{
    const bool use_parallel_assembly = (config.simulation.assembly_mode == "parallel");

    std::unique_ptr<fem::fe::ScalarFeContext> serial_context;
    std::unique_ptr<fem::fe::ParallelScalarFeContext> parallel_context;
    mfem::Mesh *mesh = nullptr;
    mfem::FiniteElementSpace *space = nullptr;

    if (use_parallel_assembly)
    {
        parallel_context = std::make_unique<fem::fe::ParallelScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh = &parallel_context->Mesh();
        space = &parallel_context->Space();
    }
    else
    {
        serial_context = std::make_unique<fem::fe::ScalarFeContext>(
            config.simulation.mesh_path, config.simulation.order,
            config.simulation.uniform_refinement_levels);
        mesh = &serial_context->Mesh();
        space = &serial_context->Space();
    }

    detail::RequireDomainMaterialsCoverAllDomains(*mesh, config, "scalar field");

    std::unordered_map<int, std::string> diffusion_by_domain;
    if (!config.domain_materials.empty())
    {
        diffusion_by_domain = detail::BuildDiffusionByDomainFromMaterials(config, materials);
    }

    const std::string diffusion_default = detail::ResolveDiffusionFromMaterial(config, materials);

    auto diffusion = detail::BuildPiecewiseDomainCoefficient(
        *mesh, diffusion_default, diffusion_by_domain, coupled_variable_names);

    auto source = detail::BuildPiecewiseDomainCoefficient(
        *mesh, field.source_default, field.source_by_domain, coupled_variable_names);

    std::unique_ptr<mfem::GridFunction> serial_solution;
    std::unique_ptr<mfem::ParGridFunction> parallel_solution;
    if (use_parallel_assembly)
    {
        auto *par_space = dynamic_cast<mfem::ParFiniteElementSpace *>(space);
        if (!par_space)
        {
            throw std::runtime_error("parallel assembly requested but ParFiniteElementSpace missing.");
        }
        parallel_solution = std::make_unique<mfem::ParGridFunction>(par_space);
        *parallel_solution = 0.0;
        detail::SolveScalarFieldOnParallelContext(*parallel_context, field, *diffusion, *source,
                                                  500, config.simulation.solver,
                                                  config.simulation.assembly_mode,
                                                  *parallel_solution);
    }
    else
    {
        serial_solution = std::make_unique<mfem::GridFunction>(space);
        *serial_solution = 0.0;
        detail::SolveScalarFieldOnContext(*serial_context, field, *diffusion, *source, 500,
                                          config.simulation.solver,
                                          config.simulation.assembly_mode, *serial_solution);
    }
    mfem::GridFunction &solution =
        use_parallel_assembly ? static_cast<mfem::GridFunction &>(*parallel_solution)
                              : static_cast<mfem::GridFunction &>(*serial_solution);

    const auto stats = detail::ComputeFieldStatistics(solution);
    spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}", field.name,
                 stats.min_value, stats.max_value, stats.mean_value, stats.nan_count);

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        ParallelAwareOutputPath(config.simulation.output_dir + "/" + field.name +
                                "_solution.txt"),
        field.name,
        *mesh, solution, "solution");

    fem::io::ParaviewExporter exporter(config.simulation.output_dir);
    exporter.Save(field.name, "solution", *mesh, solution, 0, 0.0);

    return 0;
}
}  // namespace fem::app
