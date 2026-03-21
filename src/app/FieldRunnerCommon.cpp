#include "fem/app/FieldRunnerCommon.hpp"

#include "fem/assembly/MfemLinearElasticityAssembler.hpp"
#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/fe/MechanicalFeContext.hpp"
#include "fem/fe/ScalarFeContext.hpp"
#include "fem/frontend/Expression.hpp"
#include "fem/mapping/CoefficientFactory.hpp"
#include "fem/post/MechanicalPostProcessor.hpp"
#include "fem/post/SolutionTextExporter.hpp"
#include "fem/solver/LinearSolverFactory.hpp"
#include "fem/solver/MfemPcgSolver.hpp"


#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fem::app::detail
{
namespace
{
mfem::Array<int> BuildRemainingBoundaryMarker(const mfem::Array<int> &essential_marker)
{
    mfem::Array<int> marker(essential_marker.Size());
    for (int i = 0; i < essential_marker.Size(); ++i)
    {
        marker[i] = essential_marker[i] ? 0 : 1;
    }
    return marker;
}

void MergeBoundaryMarker(mfem::Array<int> &target, const mfem::Array<int> &source)
{
    if (source.Size() == 0)
    {
        return;
    }
    if (target.Size() == 0)
    {
        target.SetSize(source.Size());
        target = 0;
    }
    const int n = std::min(target.Size(), source.Size());
    for (int i = 0; i < n; ++i)
    {
        if (source[i] != 0)
        {
            target[i] = 1;
        }
    }
}

double EvaluateAsConstant(const fem::frontend::ScalarExpression &expression)
{
    fem::frontend::EvalContext ctx;
    return expression.Evaluate(ctx);
}

std::string SelectDefaultMaterialByDomain(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials,
    const std::vector<std::string> &required_properties,
    const std::vector<std::string> &alternative_properties = {})
{
    int best_attr = std::numeric_limits<int>::max();
    std::string best_material;

    for (const auto &[material_name, domains] : config.domain_materials)
    {
        bool valid = true;
        for (const auto &property : required_properties)
        {
            if (!materials.Contains(material_name, property))
            {
                valid = false;
                break;
            }
        }

        if (!valid && !alternative_properties.empty())
        {
            valid = false;
            for (const auto &property : alternative_properties)
            {
                if (materials.Contains(material_name, property))
                {
                    valid = true;
                    break;
                }
            }
        }

        if (!valid)
        {
            continue;
        }

        for (const int attr : domains)
        {
            if (attr > 0 && attr < best_attr)
            {
                best_attr = attr;
                best_material = material_name;
            }
        }
    }

    return best_material;
}

std::vector<int> CollectMissingDomainAttributes(const mfem::Mesh &mesh,
                                                const fem::frontend::ProjectConfig &config)
{
    const int max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    std::vector<int> covered(max_attr + 1, 0);

    for (const auto &[_, domains] : config.domain_materials)
    {
        for (const int attr : domains)
        {
            if (attr > 0 && attr <= max_attr)
            {
                covered[attr] = 1;
            }
        }
    }

    std::vector<int> missing;
    for (int i = 0; i < mesh.attributes.Size(); ++i)
    {
        const int attr = mesh.attributes[i];
        if (attr > 0 && attr <= max_attr && covered[attr] == 0)
        {
            missing.push_back(attr);
        }
    }

    std::sort(missing.begin(), missing.end());
    missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
    return missing;
}

std::string FormatIntList(const std::vector<int> &values)
{
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            oss << ",";
        }
        oss << values[i];
    }
    return oss.str();
}
}  // namespace

void RequireDomainMaterialsCoverAllDomains(const mfem::Mesh &mesh,
                                           const fem::frontend::ProjectConfig &config,
                                           const std::string &context)
{
    const auto missing = CollectMissingDomainAttributes(mesh, config);
    if (!missing.empty())
    {
        throw std::runtime_error(
            "domain_materials does not fully cover mesh domain attributes for " + context +
            ". missing attributes: [" + FormatIntList(missing) + "]");
    }
}

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

ScalarBoundarySetup BuildScalarBoundarySetup(const mfem::Mesh &mesh,
                                             const fem::frontend::FieldConfig &field)
{
    ScalarBoundarySetup setup;
    setup.essential_marker = BuildBoundaryMarker(mesh, field.dirichlet_bdr_attributes, true);

    for (const auto &bc : field.dirichlet_boundary_conditions)
    {
        mfem::Array<int> marker = BuildBoundaryMarker(mesh, bc.bdr_attributes, false);
        if (marker.Size() == 0)
        {
            continue;
        }
        MergeBoundaryMarker(setup.essential_marker, marker);
        setup.dirichlet_conditions.emplace_back(bc.value, marker);
    }

    if (field.robin_on_remaining_boundaries && (field.robin_l != 0.0 || field.robin_q != 0.0))
    {
        setup.robin_conditions.emplace_back(field.robin_l, field.robin_q,
                                            BuildRemainingBoundaryMarker(setup.essential_marker));
    }
    if (!field.robin_bdr_attributes.empty() && (field.robin_l != 0.0 || field.robin_q != 0.0))
    {
        setup.robin_conditions.emplace_back(
            field.robin_l, field.robin_q,
            BuildBoundaryMarker(mesh, field.robin_bdr_attributes, false));
    }
    for (const auto &robin_bc : field.robin_boundary_conditions)
    {
        setup.robin_conditions.emplace_back(
            robin_bc.l, robin_bc.q, BuildBoundaryMarker(mesh, robin_bc.bdr_attributes, false));
    }

    return setup;
}

void SolveScalarPoissonSystem(mfem::FiniteElementSpace &space, mfem::Coefficient &diffusion,
                              mfem::Coefficient &source, const ScalarBoundarySetup &boundaries,
                              double dirichlet_value, int max_iterations,
                              const std::string &solver_name, mfem::GridFunction &solution)
{
    fem::assembly::PoissonAssemblyInput input{
        space, diffusion, source, boundaries.essential_marker, {}, {}, dirichlet_value};
    input.dirichlet_conditions = boundaries.dirichlet_conditions;
    input.robin_conditions = boundaries.robin_conditions;

    fem::assembly::MfemPoissonAssembler assembler;
    auto assembled = assembler.Assemble(input, solution);

    auto linear_solver =
        fem::solver::CreateLinearSolver(solver_name, 1e-10, 0.0, max_iterations, 0);
    linear_solver->Solve(*assembled.A.Ptr(), assembled.B, assembled.X);
    assembled.bilinear->RecoverFEMSolution(assembled.X, *assembled.linear, solution);
}

void SolveScalarFieldOnContext(fem::fe::ScalarFeContext &context,
                               const fem::frontend::FieldConfig &field,
                               mfem::Coefficient &diffusion, mfem::Coefficient &source,
                               int max_iterations, const std::string &solver_name,
                               mfem::GridFunction &solution)
{
    const ScalarBoundarySetup boundaries = BuildScalarBoundarySetup(context.Mesh(), field);
    SolveScalarPoissonSystem(context.Space(), diffusion, source, boundaries, field.dirichlet_value,
                             max_iterations, solver_name, solution);
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

std::unordered_map<int, std::string> BuildPropertyByDomainFromMaterials(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials,
    const std::string &property)
{
    std::unordered_map<int, std::string> property_by_domain;

    for (const auto &[material_name, domains] : config.domain_materials)
    {
        if (!materials.Contains(material_name, property))
        {
            continue;
        }
        const auto &expr = materials.GetProperty(material_name, property);
        const std::string raw = expr.Raw();
        for (const int domain_attr : domains)
        {
            property_by_domain[domain_attr] = raw;
        }
    }

    return property_by_domain;
}

std::string ResolveDiffusionFromMaterial(const fem::frontend::ProjectConfig &config,
                                         const fem::material::MaterialDatabase &materials)
{
    const std::string domain_default_material =
        SelectDefaultMaterialByDomain(config, materials, {"diffusion"});
    if (!domain_default_material.empty())
    {
        return materials.GetProperty(domain_default_material, "diffusion").Raw();
    }

    throw std::runtime_error("domain_materials cannot provide default diffusion material (missing "
                             "'diffusion' property).");
}

std::string ResolveElectricalConductivityFromMaterial(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials)
{
    const std::string domain_default_material = SelectDefaultMaterialByDomain(
        config, materials, {"electrical_conductivity"}, {"diffusion"});
    if (!domain_default_material.empty())
    {
        if (materials.Contains(domain_default_material, "electrical_conductivity"))
        {
            return materials.GetProperty(domain_default_material, "electrical_conductivity").Raw();
        }
        return materials.GetProperty(domain_default_material, "diffusion").Raw();
    }

    throw std::runtime_error(
        "domain_materials cannot provide default electrical conductivity material "
        "(missing 'electrical_conductivity' and fallback 'diffusion').");
}

mfem::Vector BuildPiecewiseDomainValues(const mfem::Mesh &mesh, const std::string &default_expr,
                                        const std::unordered_map<int, std::string> &overrides,
                                        const std::vector<std::string> &coupled_variable_names)
{
    const int max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    mfem::Vector values(max_attr > 0 ? max_attr : 1);

    const fem::frontend::ScalarExpression default_expression(default_expr, coupled_variable_names);
    values = EvaluateAsConstant(default_expression);

    for (const auto &[attr, expr_text] : overrides)
    {
        if (attr <= 0 || attr > max_attr)
        {
            continue;
        }
        const fem::frontend::ScalarExpression expression(expr_text, coupled_variable_names);
        values(attr - 1) = EvaluateAsConstant(expression);
    }

    return values;
}

mfem::Vector ToSizedVector(const std::vector<double> &values, int dim)
{
    mfem::Vector out(dim);
    out = 0.0;
    for (int i = 0; i < dim && i < static_cast<int>(values.size()); ++i)
    {
        out(i) = values[i];
    }
    return out;
}

LameCoefficientStorage BuildLameCoefficients(const mfem::Mesh &mesh,
                                             const fem::frontend::ProjectConfig &config,
                                             const fem::material::MaterialDatabase &materials,
                                             const std::vector<std::string> &coupled_names,
                                             const std::string &default_material)
{
    const int max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    LameCoefficientStorage storage;
    storage.lambda_values.SetSize(max_attr > 0 ? max_attr : 1);
    storage.mu_values.SetSize(max_attr > 0 ? max_attr : 1);

    auto eval_const = [&](const std::string &expr)
    {
        return fem::frontend::ScalarExpression(expr, coupled_names)
            .Evaluate(fem::frontend::EvalContext{});
    };

    double default_E = 2.1e11;
    double default_nu = 0.3;
    if (!default_material.empty())
    {
        default_E = eval_const(materials.GetProperty(default_material, "young_modulus").Raw());
        default_nu = eval_const(materials.GetProperty(default_material, "poisson_ratio").Raw());
    }

    const double default_mu = default_E / (2.0 * (1.0 + default_nu));
    const double default_lambda =
        (default_E * default_nu) / ((1.0 + default_nu) * (1.0 - 2.0 * default_nu));

    storage.lambda_values = default_lambda;
    storage.mu_values = default_mu;

    for (const auto &[material_name, domains] : config.domain_materials)
    {
        const double E = eval_const(materials.GetProperty(material_name, "young_modulus").Raw());
        const double nu = eval_const(materials.GetProperty(material_name, "poisson_ratio").Raw());
        const double mu = E / (2.0 * (1.0 + nu));
        const double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));

        for (const int attr : domains)
        {
            if (attr <= 0 || attr > storage.lambda_values.Size())
            {
                continue;
            }
            storage.lambda_values(attr - 1) = lambda;
            storage.mu_values(attr - 1) = mu;
        }
    }

    storage.lambda = std::make_unique<mfem::PWConstCoefficient>(storage.lambda_values);
    storage.mu = std::make_unique<mfem::PWConstCoefficient>(storage.mu_values);
    return storage;
}

std::unique_ptr<mfem::PWConstCoefficient> BuildAlphaCoefficientByDomain(
    const mfem::Mesh &mesh, const fem::frontend::ProjectConfig &config,
    const fem::material::MaterialDatabase &materials,
    const std::vector<std::string> &coupled_variable_names, const std::string &default_material)
{
    const int max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    mfem::Vector alpha_values(max_attr > 0 ? max_attr : 1);
    alpha_values = 0.0;

    auto eval_const = [&](const std::string &expr)
    {
        return fem::frontend::ScalarExpression(expr, coupled_variable_names)
            .Evaluate(fem::frontend::EvalContext{});
    };

    auto get_alpha_property_raw = [&](const std::string &material_name) -> std::string
    {
        if (materials.Contains(material_name, "thermal_expansion_secant"))
        {
            return materials.GetProperty(material_name, "thermal_expansion_secant").Raw();
        }
        if (materials.Contains(material_name, "thermal_expansion"))
        {
            spdlog::warn(
                "Material '{}' uses deprecated property 'thermal_expansion'; please migrate to "
                "'thermal_expansion_secant'.",
                material_name);
            return materials.GetProperty(material_name, "thermal_expansion").Raw();
        }
        throw std::runtime_error("Material '" + material_name +
                                 "' missing thermal expansion property: "
                                 "expected 'thermal_expansion_secant'.");
    };

    if (!default_material.empty())
    {
        alpha_values = eval_const(get_alpha_property_raw(default_material));
    }

    for (const auto &[material_name, domains] : config.domain_materials)
    {
        const double alpha = eval_const(get_alpha_property_raw(material_name));
        for (const int attr : domains)
        {
            if (attr <= 0 || attr > alpha_values.Size())
            {
                continue;
            }
            alpha_values(attr - 1) = alpha;
        }
    }

    return std::make_unique<mfem::PWConstCoefficient>(alpha_values);
}

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

PiecewiseVectorCoefficientStorage BuildPiecewiseDomainVectorCoefficient(
    const mfem::Mesh &mesh, int dim, const std::vector<double> &default_value,
    const std::unordered_map<int, std::vector<double>> &overrides)
{
    PiecewiseVectorCoefficientStorage storage;
    storage.coefficient = std::make_unique<mfem::PWVectorCoefficient>(dim);

    const int max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    const mfem::Vector default_vector = ToSizedVector(default_value, dim);
    for (int attr = 1; attr <= max_attr; ++attr)
    {
        auto it = overrides.find(attr);
        mfem::Vector value =
            (it != overrides.end()) ? ToSizedVector(it->second, dim) : default_vector;
        storage.pieces.push_back(std::make_unique<mfem::VectorConstantCoefficient>(value));
        storage.coefficient->UpdateCoefficient(attr, *storage.pieces.back());
    }
    return storage;
}

int RunMechanicalFieldInternal(const fem::frontend::ProjectConfig &config,
                               const fem::frontend::FieldConfig &field,
                               const std::vector<std::string> &coupled_variable_names,
                               fem::material::MaterialDatabase &materials,
                               mfem::Coefficient *temperature_by_domain)
{
    fem::fe::MechanicalFeContext fe_context(config.simulation.mesh_path, config.simulation.order,
                                            config.simulation.uniform_refinement_levels);
    mfem::Mesh &mesh = fe_context.Mesh();
    mfem::FiniteElementSpace &space = fe_context.Space();
    const int dim = fe_context.Dimension();

    RequireDomainMaterialsCoverAllDomains(mesh, config, "mechanical field");

    const auto essential_marker = BuildBoundaryMarker(mesh, field.dirichlet_bdr_attributes, true);

    std::string default_mechanical_material =
        SelectDefaultMaterialByDomain(config, materials, {"young_modulus", "poisson_ratio"});
    if (default_mechanical_material.empty())
    {
        throw std::runtime_error("domain_materials cannot provide default mechanical material "
                                 "(missing 'young_modulus'/'poisson_ratio').");
    }

    const auto lame = BuildLameCoefficients(mesh, config, materials, coupled_variable_names,
                                            default_mechanical_material);
    auto body_force = BuildPiecewiseDomainVectorCoefficient(mesh, dim, field.body_force_default,
                                                            field.body_force_by_domain);

    const mfem::Vector dirichlet_value = ToSizedVector(field.dirichlet_displacement_value, dim);
    mfem::VectorConstantCoefficient dirichlet_disp(dirichlet_value);

    std::vector<fem::assembly::MechanicalTractionBoundary> traction_boundaries;
    for (const auto &traction_bc : field.traction_boundary_conditions)
    {
        mfem::Array<int> marker = BuildBoundaryMarker(mesh, traction_bc.bdr_attributes, false);
        if (marker.Size() == 0)
        {
            continue;
        }
        traction_boundaries.emplace_back(ToSizedVector(traction_bc.value, dim), marker);
    }

    std::vector<fem::assembly::MechanicalPressureBoundary> pressure_boundaries;
    for (const auto &pressure_bc : field.pressure_boundary_conditions)
    {
        mfem::Array<int> marker = BuildBoundaryMarker(mesh, pressure_bc.bdr_attributes, false);
        if (marker.Size() == 0)
        {
            continue;
        }
        pressure_boundaries.emplace_back(-pressure_bc.value, marker);
    }

    std::vector<fem::assembly::MechanicalNormalDisplacementBoundary> normal_displacement_boundaries;
    for (const auto &normal_bc : field.normal_displacement_boundary_conditions)
    {
        mfem::Array<int> marker = BuildBoundaryMarker(mesh, normal_bc.bdr_attributes, false);
        if (marker.Size() == 0)
        {
            continue;
        }
        normal_displacement_boundaries.emplace_back(normal_bc.penalty, marker);
    }

    std::unique_ptr<mfem::PWConstCoefficient> alpha_by_domain;
    std::unique_ptr<mfem::H1_FECollection> thermal_interp_fec;
    std::unique_ptr<mfem::FiniteElementSpace> thermal_interp_space;
    std::unique_ptr<mfem::GridFunction> thermal_interp;
    std::unique_ptr<mfem::GridFunctionCoefficient> thermal_interp_coeff;
    mfem::Coefficient *temperature_for_mechanical = nullptr;
    if (temperature_by_domain && field.enable_thermal_strain_coupling)
    {
        std::string default_thermal_expansion_material = SelectDefaultMaterialByDomain(
            config, materials, {"thermal_expansion_secant"}, {"thermal_expansion"});
        if (default_thermal_expansion_material.empty())
        {
            throw std::runtime_error(
                "domain_materials cannot provide default thermal expansion material "
                "(missing 'thermal_expansion_secant' and fallback 'thermal_expansion').");
        }

        alpha_by_domain = BuildAlphaCoefficientByDomain(
            mesh, config, materials, coupled_variable_names, default_thermal_expansion_material);

        thermal_interp_fec = std::make_unique<mfem::H1_FECollection>(config.simulation.order, dim);
        thermal_interp_space =
            std::make_unique<mfem::FiniteElementSpace>(&mesh, thermal_interp_fec.get());
        thermal_interp = std::make_unique<mfem::GridFunction>(thermal_interp_space.get());
        thermal_interp->ProjectCoefficient(*temperature_by_domain);
        thermal_interp_coeff =
            std::make_unique<mfem::GridFunctionCoefficient>(thermal_interp.get());
        temperature_for_mechanical = thermal_interp_coeff.get();
    }

    mfem::GridFunction displacement(&space);
    displacement = 0.0;

    fem::assembly::LinearElasticityInput input{space,
                                               *lame.lambda,
                                               *lame.mu,
                                               *body_force.coefficient,
                                               essential_marker,
                                               dirichlet_disp,
                                               traction_boundaries,
                                               {},
                                               {},
                                               nullptr,
                                               nullptr,
                                               field.reference_temperature};
    input.pressure_boundaries = std::move(pressure_boundaries);
    input.normal_displacement_boundaries = std::move(normal_displacement_boundaries);
    input.reference_temperature = field.reference_temperature;
    if (alpha_by_domain && temperature_for_mechanical)
    {
        input.thermal_expansion_secant = alpha_by_domain.get();
        input.temperature = temperature_for_mechanical;
    }

    fem::assembly::MfemLinearElasticityAssembler assembler;
    auto assembled = assembler.Assemble(input, displacement);

    auto linear_solver =
        fem::solver::CreateLinearSolver(config.simulation.solver, 1e-10, 0.0, 2000, 0);
    linear_solver->Solve(*assembled.A.Ptr(), assembled.B, assembled.X);
    assembled.bilinear->RecoverFEMSolution(assembled.X, *assembled.linear, displacement);

    mfem::H1_FECollection von_fec(1, dim);
    mfem::FiniteElementSpace von_space(&mesh, &von_fec);
    mfem::GridFunction von_mises(&von_space);

    std::unique_ptr<fem::post::ThermalStrainParameters> thermal_post;
    if (alpha_by_domain && temperature_for_mechanical)
    {
        thermal_post =
            std::make_unique<fem::post::ThermalStrainParameters>(fem::post::ThermalStrainParameters{
                *alpha_by_domain, *temperature_for_mechanical, field.reference_temperature});
    }

    if (field.enable_stress_postprocess)
    {
        fem::post::MechanicalPostProcessor::FillVonMisesNodalField(
            mesh, displacement, fem::post::LameParameters{*lame.lambda, *lame.mu}, von_space,
            von_mises, thermal_post.get());

        const std::string csv_path =
            config.simulation.output_dir + "/" + field.name + "_element_stress.csv";
        fem::post::MechanicalPostProcessor::ExportElementStressCsv(
            csv_path, mesh, displacement, fem::post::LameParameters{*lame.lambda, *lame.mu},
            thermal_post.get());

        fem::post::SolutionTextExporter::ExportScalarNodalTxt(
            config.simulation.output_dir + "/" + field.name + "_von_mises.txt", field.name, mesh,
            von_mises, "von_mises");
    }

    fem::post::SolutionTextExporter::ExportVectorNodalTxt(
        config.simulation.output_dir + "/" + field.name + "_solution.txt", field.name, mesh,
        displacement, "u");

    std::filesystem::create_directories(config.simulation.output_dir);
    mfem::ParaViewDataCollection collection(field.name, &mesh);
    collection.SetPrefixPath(config.simulation.output_dir);
    collection.SetLevelsOfDetail(1);
    collection.SetCycle(0);
    collection.SetTime(0.0);
    collection.RegisterField("displacement", &displacement);
    if (field.enable_stress_postprocess)
    {
        collection.RegisterField("von_mises", &von_mises);
    }
    collection.Save();

    return 0;
}
}  // namespace fem::app::detail
