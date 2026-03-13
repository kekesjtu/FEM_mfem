#include "fem/app/Application.hpp"

#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/fe/FeContext.hpp"
#include "fem/frontend/ConfigLoader.hpp"
#include "fem/frontend/Expression.hpp"
#include "fem/io/ParaviewExporter.hpp"
#include "fem/log/Logger.hpp"
#include "fem/mapping/CoefficientFactory.hpp"
#include "fem/material/MaterialDatabase.hpp"
#include "fem/solver/MfemPcgSolver.hpp"

#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace fem::app
{
namespace
{
mfem::Array<int> BuildBoundaryMarker(const mfem::Mesh &mesh, const std::vector<int> &attributes,
                                     bool fallback_to_all_if_invalid_selection)
{
    mfem::Array<int> marker;
    int max_bdr_attr = 0;
    if (mesh.bdr_attributes.Size() > 0)
    {
        max_bdr_attr = mesh.bdr_attributes.Max();
    }
    else if (mesh.GetNBE() > 0)
    {
        max_bdr_attr = 1;
    }
    marker.SetSize(max_bdr_attr);
    marker = 0;
    int marked_count = 0;

    for (const int attr : attributes)
    {
        if (attr <= 0 || attr > max_bdr_attr)
        {
            continue;
        }
        if (marker[attr - 1] == 0)
        {
            marker[attr - 1] = 1;
            ++marked_count;
        }
    }

    if (max_bdr_attr > 0 && !attributes.empty() && marked_count == 0 &&
        fallback_to_all_if_invalid_selection)
    {
        marker = 1;
        spdlog::warn(
            "No valid Dirichlet boundary attribute provided; fallback to all boundary attributes.");
    }

    return marker;
}

mfem::Array<int> BuildRemainingBoundaryMarker(const mfem::Array<int> &essential_marker)
{
    mfem::Array<int> marker(essential_marker.Size());
    for (int i = 0; i < essential_marker.Size(); ++i)
    {
        marker[i] = essential_marker[i] ? 0 : 1;
    }
    return marker;
}

double EvaluateAsConstant(const fem::frontend::ScalarExpression &expression)
{
    fem::frontend::EvalContext ctx;
    return expression.Evaluate(ctx);
}

std::unique_ptr<mfem::Coefficient> BuildPiecewiseDomainCoefficient(
    const mfem::Mesh &mesh, const std::string &default_expr,
    const std::unordered_map<int, std::string> &overrides,
    const std::vector<std::string> &coupled_variable_names)
{
    if (mesh.attributes.Size() == 0 || overrides.empty())
    {
        return fem::mapping::CoefficientFactory::CreateScalar(
            fem::frontend::ScalarExpression(default_expr, coupled_variable_names));
    }

    const int max_attr = mesh.attributes.Max();
    mfem::Vector attr_values(max_attr);

    const fem::frontend::ScalarExpression default_expression(default_expr, coupled_variable_names);
    attr_values = EvaluateAsConstant(default_expression);

    for (const auto &[attr, expr_text] : overrides)
    {
        if (attr <= 0 || attr > max_attr)
        {
            continue;
        }
        const fem::frontend::ScalarExpression expression(expr_text, coupled_variable_names);
        attr_values(attr - 1) = EvaluateAsConstant(expression);
    }

    return std::make_unique<mfem::PWConstCoefficient>(attr_values);
}

std::unordered_map<int, std::string> BuildDiffusionByDomainFromMaterials(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials)
{
    std::unordered_map<int, std::string> diffusion_by_domain;

    for (const auto &[material_name, domains] : config.domain_materials)
    {
        const auto &expr = materials.GetProperty(material_name, "diffusion");
        const std::string raw = expr.Raw();
        for (const int domain_attr : domains)
        {
            diffusion_by_domain[domain_attr] = raw;
        }
    }

    return diffusion_by_domain;
}

struct FieldStatistics
{
    double min_value = 0.0;
    double max_value = 0.0;
    double mean_value = 0.0;
    int nan_count = 0;
};

FieldStatistics ComputeFieldStatistics(const mfem::GridFunction &field)
{
    FieldStatistics stats;
    if (field.Size() == 0)
    {
        return stats;
    }

    stats.min_value = field(0);
    stats.max_value = field(0);
    double sum = 0.0;

    for (int i = 0; i < field.Size(); ++i)
    {
        const double value = field(i);
        if (!std::isfinite(value))
        {
            ++stats.nan_count;
            continue;
        }
        if (value < stats.min_value)
        {
            stats.min_value = value;
        }
        if (value > stats.max_value)
        {
            stats.max_value = value;
        }
        sum += value;
    }

    stats.mean_value = sum / static_cast<double>(field.Size() - stats.nan_count);
    return stats;
}
}  // namespace

int Application::Run(const std::string &config_path)
{
    const auto config = fem::frontend::ConfigLoader::LoadFromFile(config_path);
    fem::log::Init(config.simulation.log_level);

    spdlog::info("Loaded configuration: {}", config_path);
    spdlog::info("Mesh: {}, order: {}, refinement: {}", config.simulation.mesh_path,
                 config.simulation.order, config.simulation.uniform_refinement_levels);

    if (config.fields.empty())
    {
        throw std::runtime_error("No field definitions in config.");
    }

    std::vector<std::string> coupled_variable_names;
    coupled_variable_names.reserve(config.fields.size());
    for (const auto &configured_field : config.fields)
    {
        coupled_variable_names.push_back(configured_field.name);
    }

    const auto &field = config.fields.front();
    const bool supported_type = (field.type == "electrostatic") || (field.type == "thermal") ||
                                (field.type == "mechanical");
    if (!supported_type)
    {
        throw std::runtime_error("Unsupported field type: " + field.type);
    }

    fem::material::MaterialDatabase materials;
    for (const auto &material : config.materials)
    {
        materials.Add(material, coupled_variable_names);
    }

    fem::fe::FeContext fe_context(config.simulation.mesh_path, config.simulation.order,
                                  config.simulation.uniform_refinement_levels);

    spdlog::info("Running field '{}' of type '{}'.", field.name, field.type);

    std::unordered_map<int, std::string> diffusion_by_domain = field.diffusion_by_domain;
    if (diffusion_by_domain.empty() && !config.domain_materials.empty())
    {
        diffusion_by_domain = BuildDiffusionByDomainFromMaterials(config, materials);
    }

    std::string diffusion_default = field.diffusion;
    if (diffusion_default == "1.0" && !field.material.empty())
    {
        diffusion_default = materials.GetProperty(field.material, "diffusion").Raw();
    }

    auto diffusion = BuildPiecewiseDomainCoefficient(fe_context.Mesh(), diffusion_default,
                                                     diffusion_by_domain, coupled_variable_names);

    auto source = BuildPiecewiseDomainCoefficient(fe_context.Mesh(), field.source_default,
                                                  field.source_by_domain, coupled_variable_names);

    const auto essential_marker =
        BuildBoundaryMarker(fe_context.Mesh(), field.dirichlet_bdr_attributes, true);

    mfem::Array<int> robin_marker;
    robin_marker.SetSize(essential_marker.Size());
    robin_marker = 0;
    if (field.robin_on_remaining_boundaries)
    {
        robin_marker = BuildRemainingBoundaryMarker(essential_marker);
    }
    if (!field.robin_bdr_attributes.empty())
    {
        robin_marker = BuildBoundaryMarker(fe_context.Mesh(), field.robin_bdr_attributes, false);
    }

    std::unique_ptr<mfem::ConstantCoefficient> robin_l;
    std::unique_ptr<mfem::ConstantCoefficient> robin_q;
    if (field.robin_l != 0.0 || field.robin_q != 0.0)
    {
        robin_l = std::make_unique<mfem::ConstantCoefficient>(field.robin_l);
        robin_q = std::make_unique<mfem::ConstantCoefficient>(field.robin_q);
    }

    mfem::GridFunction solution(&fe_context.Space());
    solution = 0.0;

    fem::assembly::PoissonAssemblyInput input{
        fe_context.Space(), *diffusion,    *source,      essential_marker,
        robin_l.get(),      robin_q.get(), robin_marker, field.dirichlet_value};

    fem::assembly::MfemPoissonAssembler assembler;
    auto assembled = assembler.Assemble(input, solution);

    fem::solver::MfemPcgSolver linear_solver(1e-10, 0.0, 500, 0);
    linear_solver.Solve(*assembled.A.Ptr(), assembled.B, assembled.X);
    assembled.bilinear->RecoverFEMSolution(assembled.X, *assembled.linear, solution);

    const auto stats = ComputeFieldStatistics(solution);
    spdlog::info("Field '{}' stats: min={}, max={}, mean={}, non_finite_count={}", field.name,
                 stats.min_value, stats.max_value, stats.mean_value, stats.nan_count);

    fem::io::ParaviewExporter exporter(config.simulation.output_dir);
    exporter.Save(field.name, "solution", fe_context.Mesh(), solution, 0, 0.0);

    spdlog::info("Finished field '{}' and wrote ParaView output to '{}'.", field.name,
                 config.simulation.output_dir);
    return 0;
}
}  // namespace fem::app
