#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/material/MaterialDatabase.hpp"

#include "fem/assembly/IPoissonAssembler.hpp"

#include "mfem.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem::fe
{
class ScalarFeContext;
}  // namespace fem::fe

namespace fem::app::detail
{
struct ScalarBoundarySetup
{
    mfem::Array<int> essential_marker;
    std::vector<fem::assembly::DirichletBoundaryCondition> dirichlet_conditions;
    std::vector<fem::assembly::RobinBoundaryCondition> robin_conditions;
};

struct LameCoefficientStorage
{
    mfem::Vector lambda_values;
    mfem::Vector mu_values;
    std::unique_ptr<mfem::PWConstCoefficient> lambda;
    std::unique_ptr<mfem::PWConstCoefficient> mu;
};

struct PiecewiseVectorCoefficientStorage
{
    std::vector<std::unique_ptr<mfem::VectorConstantCoefficient>> pieces;
    std::unique_ptr<mfem::PWVectorCoefficient> coefficient;
};

struct FieldStatistics
{
    double min_value = 0.0;
    double max_value = 0.0;
    double mean_value = 0.0;
    int nan_count = 0;
};

mfem::Array<int> BuildBoundaryMarker(const mfem::Mesh &mesh, const std::vector<int> &attributes,
                                     bool fallback_to_all_if_invalid_selection);

ScalarBoundarySetup BuildScalarBoundarySetup(const mfem::Mesh &mesh,
                                             const fem::frontend::FieldConfig &field);

void SolveScalarPoissonSystem(mfem::FiniteElementSpace &space, mfem::Coefficient &diffusion,
                              mfem::Coefficient &source, const ScalarBoundarySetup &boundaries,
                              double dirichlet_value, int max_iterations,
                              const std::string &solver_name, mfem::GridFunction &solution);

void SolveScalarFieldOnContext(fem::fe::ScalarFeContext &context,
                               const fem::frontend::FieldConfig &field,
                               mfem::Coefficient &diffusion, mfem::Coefficient &source,
                               int max_iterations, const std::string &solver_name,
                               mfem::GridFunction &solution);

std::unique_ptr<mfem::Coefficient> BuildPiecewiseDomainCoefficient(
    const mfem::Mesh &mesh, const std::string &default_expr,
    const std::unordered_map<int, std::string> &overrides,
    const std::vector<std::string> &coupled_variable_names);

std::unordered_map<int, std::string> BuildDiffusionByDomainFromMaterials(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials);

std::unordered_map<int, std::string> BuildPropertyByDomainFromMaterials(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials,
    const std::string &property);

void RequireDomainMaterialsCoverAllDomains(const mfem::Mesh &mesh,
                                           const fem::frontend::ProjectConfig &config,
                                           const std::string &context);

std::string ResolveDiffusionFromMaterial(const fem::frontend::ProjectConfig &config,
                                         const fem::material::MaterialDatabase &materials);

std::string ResolveElectricalConductivityFromMaterial(
    const fem::frontend::ProjectConfig &config, const fem::material::MaterialDatabase &materials);

mfem::Vector BuildPiecewiseDomainValues(const mfem::Mesh &mesh, const std::string &default_expr,
                                        const std::unordered_map<int, std::string> &overrides,
                                        const std::vector<std::string> &coupled_variable_names);

mfem::Vector ToSizedVector(const std::vector<double> &values, int dim);

LameCoefficientStorage BuildLameCoefficients(const mfem::Mesh &mesh,
                                             const fem::frontend::ProjectConfig &config,
                                             const fem::material::MaterialDatabase &materials,
                                             const std::vector<std::string> &coupled_names,
                                             const std::string &default_material);

std::unique_ptr<mfem::PWConstCoefficient> BuildAlphaCoefficientByDomain(
    const mfem::Mesh &mesh, const fem::frontend::ProjectConfig &config,
    const fem::material::MaterialDatabase &materials,
    const std::vector<std::string> &coupled_variable_names, const std::string &default_material);

FieldStatistics ComputeFieldStatistics(const mfem::GridFunction &field);

PiecewiseVectorCoefficientStorage BuildPiecewiseDomainVectorCoefficient(
    const mfem::Mesh &mesh, int dim, const std::vector<double> &default_value,
    const std::unordered_map<int, std::vector<double>> &overrides);

int RunMechanicalFieldInternal(const fem::frontend::ProjectConfig &config,
                               const fem::frontend::FieldConfig &field,
                               const std::vector<std::string> &coupled_variable_names,
                               fem::material::MaterialDatabase &materials,
                               mfem::Coefficient *temperature_by_domain);
}  // namespace fem::app::detail
