#pragma once

#include "fem/frontend/Config.hpp"

#include <string>

namespace fem::frontend
{
class ConfigLoader
{
  public:
    static ProjectConfig LoadFromFile(const std::string &path);

  private:
    ProjectConfig config_;

  private:
    template <typename T>
    T ReadOrDefault(const nlohmann::json &node, const std::string &key, T fallback);

    MaterialDatabase LoadMaterialDatabase(const nlohmann::json &materials_node);
    ScalarFieldConfig LoadScalarFieldConfig(const nlohmann::json &fields_node);
    MechanicalFieldConfig LoadMechanicalFieldConfig(const nlohmann::json &fields_node);

    FEconfig BuildFEconfig();// 这个函数负责根据config_中的信息构建FEconfig结构体，包含网格、有限元集合和有限元空间的指针
};
}  // namespace fem::frontend
