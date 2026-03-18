#include "fem/frontend/ConfigLoader.hpp"

#include "nlohmann/json.hpp"

#include <fstream>
#include <stdexcept>

namespace fem::frontend
{
using json = nlohmann::json;

namespace
{
template <typename T>
T ReadOrDefault(const json &node, const std::string &key, T fallback)
{
    if (!node.contains(key))
    {
        return fallback;
    }
    return node.at(key).get<T>();
}

std::unordered_map<int, std::string> ParseIntStringMap(const json &node)
{
    std::unordered_map<int, std::string> result;
    for (const auto &[k, v] : node.items())
    {
        result.emplace(std::stoi(k), v.get<std::string>());
    }
    return result;
}

std::unordered_map<int, std::vector<double>> ParseIntVectorMap(const json &node)
{
    std::unordered_map<int, std::vector<double>> result;
    for (const auto &[k, v] : node.items())
    {
        result.emplace(std::stoi(k), v.get<std::vector<double>>());
    }
    return result;
}
}  // namespace

ProjectConfig ConfigLoader::LoadFromFile(const std::string &path)
{
    std::ifstream input(path);
    if (!input)
    {
        throw std::runtime_error("Failed to open config file: " + path);
    }

    json root;
    input >> root;

    ProjectConfig config;

    const auto &sim = root.at("simulation");
    config.simulation.mesh_path = sim.at("mesh_path").get<std::string>();
    config.simulation.order = ReadOrDefault<int>(sim, "order", 1);
    config.simulation.uniform_refinement_levels =
        ReadOrDefault<int>(sim, "uniform_refinement_levels", 0);
    config.simulation.output_dir = ReadOrDefault<std::string>(sim, "output_dir", "results");
    config.simulation.log_level = ReadOrDefault<std::string>(sim, "log_level", "info");

    if (root.contains("materials"))
    {
        for (const auto &material_node : root.at("materials"))
        {
            MaterialConfig material;
            material.name = material_node.at("name").get<std::string>();

            if (material_node.contains("properties"))
            {
                for (const auto &[property, value] : material_node.at("properties").items())
                {
                    material.properties.emplace(property, value.get<std::string>());
                }
            }

            config.materials.push_back(std::move(material));
        }
    }

    if (root.contains("domain_materials"))
    {
        for (const auto &[material_name, domains] : root.at("domain_materials").items())
        {
            config.domain_materials.emplace(material_name, domains.get<std::vector<int>>());
        }
    }

    if (root.contains("fields"))
    {
        for (const auto &field_node : root.at("fields"))
        {
            FieldConfig field;
            field.name = field_node.at("name").get<std::string>();
            field.type = field_node.at("type").get<std::string>();
            field.material = ReadOrDefault<std::string>(field_node, "material", "");
            field.source = ReadOrDefault<std::string>(field_node, "source", "0.0");
            field.source_default =
                ReadOrDefault<std::string>(field_node, "source_default", field.source);
            if (field_node.contains("source_by_domain"))
            {
                field.source_by_domain = ParseIntStringMap(field_node.at("source_by_domain"));
            }
            field.dirichlet_value = ReadOrDefault<double>(field_node, "dirichlet_value", 0.0);
            field.dirichlet_bdr_attributes = ReadOrDefault<std::vector<int>>(
                field_node, "dirichlet_bdr_attributes", std::vector<int>{});
            if (field_node.contains("dirichlet_boundary_conditions"))
            {
                for (const auto &bc_node : field_node.at("dirichlet_boundary_conditions"))
                {
                    FieldConfig::ScalarDirichletBoundaryCondition bc;
                    bc.bdr_attributes = ReadOrDefault<std::vector<int>>(bc_node, "bdr_attributes",
                                                                        std::vector<int>{});
                    bc.value = ReadOrDefault<double>(bc_node, "value", 0.0);
                    field.dirichlet_boundary_conditions.push_back(std::move(bc));
                }
            }
            field.robin_l = ReadOrDefault<double>(field_node, "robin_l", 0.0);
            field.robin_q = ReadOrDefault<double>(field_node, "robin_q", 0.0);
            field.robin_bdr_attributes = ReadOrDefault<std::vector<int>>(
                field_node, "robin_bdr_attributes", std::vector<int>{});
            field.robin_on_remaining_boundaries =
                ReadOrDefault<bool>(field_node, "robin_on_remaining_boundaries", false);
            if (field_node.contains("robin_boundary_conditions"))
            {
                for (const auto &bc_node : field_node.at("robin_boundary_conditions"))
                {
                    FieldConfig::RobinBoundaryCondition bc;
                    bc.bdr_attributes = ReadOrDefault<std::vector<int>>(bc_node, "bdr_attributes",
                                                                        std::vector<int>{});
                    bc.l = ReadOrDefault<double>(bc_node, "l", 0.0);
                    bc.q = ReadOrDefault<double>(bc_node, "q", 0.0);
                    field.robin_boundary_conditions.push_back(std::move(bc));
                }
            }

            field.body_force_default = ReadOrDefault<std::vector<double>>(
                field_node, "body_force_default", std::vector<double>{0.0, 0.0, 0.0});
            if (field_node.contains("body_force_by_domain"))
            {
                field.body_force_by_domain =
                    ParseIntVectorMap(field_node.at("body_force_by_domain"));
            }

            field.dirichlet_displacement_value = ReadOrDefault<std::vector<double>>(
                field_node, "dirichlet_displacement_value", std::vector<double>{0.0, 0.0, 0.0});

            if (field_node.contains("traction_boundary_conditions"))
            {
                for (const auto &bc_node : field_node.at("traction_boundary_conditions"))
                {
                    FieldConfig::VectorBoundaryCondition bc;
                    bc.bdr_attributes = ReadOrDefault<std::vector<int>>(bc_node, "bdr_attributes",
                                                                        std::vector<int>{});
                    bc.value = ReadOrDefault<std::vector<double>>(
                        bc_node, "value", std::vector<double>{0.0, 0.0, 0.0});
                    field.traction_boundary_conditions.push_back(std::move(bc));
                }
            }

            if (field_node.contains("pressure_boundary_conditions"))
            {
                for (const auto &bc_node : field_node.at("pressure_boundary_conditions"))
                {
                    FieldConfig::ScalarBoundaryCondition bc;
                    bc.bdr_attributes = ReadOrDefault<std::vector<int>>(bc_node, "bdr_attributes",
                                                                        std::vector<int>{});
                    bc.value = ReadOrDefault<double>(bc_node, "value", 0.0);
                    field.pressure_boundary_conditions.push_back(std::move(bc));
                }
            }

            if (field_node.contains("normal_displacement_boundary_conditions"))
            {
                for (const auto &bc_node : field_node.at("normal_displacement_boundary_conditions"))
                {
                    FieldConfig::NormalDisplacementBoundaryCondition bc;
                    bc.bdr_attributes = ReadOrDefault<std::vector<int>>(bc_node, "bdr_attributes",
                                                                        std::vector<int>{});
                    bc.penalty = ReadOrDefault<double>(bc_node, "penalty", 1.0e12);
                    field.normal_displacement_boundary_conditions.push_back(std::move(bc));
                }
            }

            field.reference_temperature =
                ReadOrDefault<double>(field_node, "reference_temperature", 293.15);
            field.enable_thermal_strain_coupling =
                ReadOrDefault<bool>(field_node, "enable_thermal_strain_coupling", true);

            field.enable_stress_postprocess =
                ReadOrDefault<bool>(field_node, "enable_stress_postprocess", true);
            config.fields.push_back(std::move(field));
        }
    }

    return config;
}
}  // namespace fem::frontend
