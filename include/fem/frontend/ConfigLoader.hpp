#pragma once

#include "fem/frontend/Config.hpp"
#include "nlohmann/json.hpp"

#include <string>

namespace fem::frontend
{
class ConfigLoader
{
  public:
    static ProjectConfig LoadFromFile(const std::string &path);

  private:
    void LoadSimulation(const nlohmann::json &root);
    void LoadMaterials(const nlohmann::json &root);
    void LoadFields(const nlohmann::json &root);
    ScalarFieldConfig LoadScalarField(const nlohmann::json &field_node);
    MechanicalFieldConfig LoadMechanicalField(const nlohmann::json &field_node);
    void BuildFE();

    ProjectConfig config_;
};
}  // namespace fem::frontend
