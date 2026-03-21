#include "fem/app/FieldRunnerCommon.hpp"

#include "fem/app/FieldRunners.hpp"
#include "fem/fe/ScalarFeContext.hpp"
#include "fem/io/ParaviewExporter.hpp"
#include "fem/post/SolutionTextExporter.hpp"

#include "spdlog/spdlog.h"

namespace fem::app
{
int RunScalarField(const fem::frontend::ProjectConfig &config,
                   const fem::frontend::FieldConfig &field,
                   const std::vector<std::string> &coupled_variable_names,
                   fem::material::MaterialDatabase &materials)
{
    fem::fe::ScalarFeContext fe_context(config.simulation.mesh_path, config.simulation.order,
                                        config.simulation.uniform_refinement_levels);

    detail::RequireDomainMaterialsCoverAllDomains(fe_context.Mesh(), config, "scalar field");

    std::unordered_map<int, std::string> diffusion_by_domain;
    if (!config.domain_materials.empty())
    {
        diffusion_by_domain = detail::BuildDiffusionByDomainFromMaterials(config, materials);
    }

    const std::string diffusion_default = detail::ResolveDiffusionFromMaterial(config, materials);

    auto diffusion = detail::BuildPiecewiseDomainCoefficient(
        fe_context.Mesh(), diffusion_default, diffusion_by_domain, coupled_variable_names);

    auto source = detail::BuildPiecewiseDomainCoefficient(
        fe_context.Mesh(), field.source_default, field.source_by_domain, coupled_variable_names);

    mfem::GridFunction solution(&fe_context.Space());
    solution = 0.0;
    detail::SolveScalarFieldOnContext(fe_context, field, *diffusion, *source, 500,
                                      config.simulation.solver, solution);

    const auto stats = detail::ComputeFieldStatistics(solution);
    spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}", field.name,
                 stats.min_value, stats.max_value, stats.mean_value, stats.nan_count);

    fem::post::SolutionTextExporter::ExportScalarNodalTxt(
        config.simulation.output_dir + "/" + field.name + "_solution.txt", field.name,
        fe_context.Mesh(), solution, "solution");

    fem::io::ParaviewExporter exporter(config.simulation.output_dir);
    exporter.Save(field.name, "solution", fe_context.Mesh(), solution, 0, 0.0);

    return 0;
}
}  // namespace fem::app
