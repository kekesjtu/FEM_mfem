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
            field.diffusion = ReadOrDefault<std::string>(field_node, "diffusion", "1.0");
            if (field_node.contains("diffusion_by_domain"))
            {
                field.diffusion_by_domain = ParseIntStringMap(field_node.at("diffusion_by_domain"));
            }
            field.dirichlet_value = ReadOrDefault<double>(field_node, "dirichlet_value", 0.0);
            field.dirichlet_bdr_attributes = ReadOrDefault<std::vector<int>>(
                field_node, "dirichlet_bdr_attributes", std::vector<int>{});
            field.robin_l = ReadOrDefault<double>(field_node, "robin_l", 0.0);
            field.robin_q = ReadOrDefault<double>(field_node, "robin_q", 0.0);
            field.robin_bdr_attributes = ReadOrDefault<std::vector<int>>(
                field_node, "robin_bdr_attributes", std::vector<int>{});
            field.robin_on_remaining_boundaries =
                ReadOrDefault<bool>(field_node, "robin_on_remaining_boundaries", false);
            config.fields.push_back(std::move(field));
        }
    }

    return config;
}
}  // namespace fem::frontend
