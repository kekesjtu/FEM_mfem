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
};

struct MaterialConfig
{
    std::string name;
    std::unordered_map<std::string, std::string> properties;
};

struct FieldConfig
{
    std::string name;
    std::string type;
    std::string material;
    std::string source = "0.0";
    std::string source_default = "0.0";
    std::unordered_map<int, std::string> source_by_domain;
    std::string diffusion = "1.0";
    std::unordered_map<int, std::string> diffusion_by_domain;
    double dirichlet_value = 0.0;
    std::vector<int> dirichlet_bdr_attributes;
    double robin_l = 0.0;
    double robin_q = 0.0;
    std::vector<int> robin_bdr_attributes;
    bool robin_on_remaining_boundaries = false;
};

struct ProjectConfig
{
    SimulationConfig simulation;
    std::vector<MaterialConfig> materials;
    std::unordered_map<std::string, std::vector<int>> domain_materials;
    std::vector<FieldConfig> fields;
};
}  // namespace fem::frontend
