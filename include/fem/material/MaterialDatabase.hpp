#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/frontend/Expression.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem::material
{
class MaterialDatabase
{
  public:
    void Add(const frontend::MaterialConfig &config,
             const std::vector<std::string> &coupled_variable_names = {})
    {
        std::unordered_map<std::string, frontend::ScalarExpression> property_map;
        property_map.reserve(config.properties.size());
        for (const auto &[key, expr] : config.properties)
        {
            property_map.emplace(key, frontend::ScalarExpression(expr, coupled_variable_names));
        }
        materials_.emplace(config.name, std::move(property_map));
    }

    const frontend::ScalarExpression &GetProperty(const std::string &material,
                                                  const std::string &property) const
    {
        const auto mat_it = materials_.find(material);
        if (mat_it == materials_.end())
        {
            throw std::runtime_error("Material not found: " + material);
        }

        const auto prop_it = mat_it->second.find(property);
        if (prop_it == mat_it->second.end())
        {
            throw std::runtime_error("Property not found: " + material + "." + property);
        }

        return prop_it->second;
    }

  private:
    std::unordered_map<std::string, std::unordered_map<std::string, frontend::ScalarExpression>>
        materials_;
};
}  // namespace fem::material
