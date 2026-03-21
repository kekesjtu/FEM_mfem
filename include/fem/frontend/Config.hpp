#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace fem::frontend
{
struct SimulationConfig
{
    std::string mesh_path;
    int order = 1;
    int uniform_refinement_levels = 0;
    std::string output_dir = "results";
    std::string log_level = "info";
    std::string solver = "pcg";
};

struct MaterialConfig
{
    std::string name;
    std::unordered_map<std::string, std::string> properties;
};

struct FieldConfig
{
    struct VectorBoundaryCondition
    {
        std::vector<int> bdr_attributes;
        std::vector<double> value;
    };

    struct ScalarBoundaryCondition
    {
        std::vector<int> bdr_attributes;
        double value = 0.0;
    };

    struct NormalDisplacementBoundaryCondition
    {
        std::vector<int> bdr_attributes;
        double penalty = 1.0e12;
    };

    struct RobinBoundaryCondition
    {
        std::vector<int> bdr_attributes;
        double l = 0.0;
        double q = 0.0;
    };

    struct ScalarDirichletBoundaryCondition
    {
        std::vector<int> bdr_attributes;
        double value = 0.0;
    };

    std::string name;
    std::string type;
    std::string material;
    std::string source = "0.0";
    std::string source_default = "0.0";
    std::unordered_map<int, std::string> source_by_domain;
    double dirichlet_value = 0.0;
    std::vector<int> dirichlet_bdr_attributes;
    std::vector<ScalarDirichletBoundaryCondition> dirichlet_boundary_conditions;
    double robin_l = 0.0;
    double robin_q = 0.0;
    std::vector<int> robin_bdr_attributes;
    bool robin_on_remaining_boundaries = false;
    std::vector<RobinBoundaryCondition> robin_boundary_conditions;

    std::vector<double> body_force_default{0.0, 0.0, 0.0};
    std::unordered_map<int, std::vector<double>> body_force_by_domain;
    std::vector<double> dirichlet_displacement_value{0.0, 0.0, 0.0};
    std::vector<NormalDisplacementBoundaryCondition> normal_displacement_boundary_conditions;
    std::vector<VectorBoundaryCondition> traction_boundary_conditions;
    std::vector<ScalarBoundaryCondition> pressure_boundary_conditions;
    double reference_temperature = 293.15;
    bool enable_thermal_strain_coupling = true;
    bool enable_stress_postprocess = true;
};

struct ProjectConfig
{
    SimulationConfig simulation;
    std::vector<MaterialConfig> materials;
    std::unordered_map<std::string, std::vector<int>> domain_materials;
    std::vector<FieldConfig> fields;
};
}  // namespace fem::frontend
