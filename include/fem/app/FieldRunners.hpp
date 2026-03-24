#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/material/MaterialDatabase.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem::app
{
using FieldRunner =
    std::function<int(const fem::frontend::ProjectConfig &, const fem::frontend::FieldConfig &,
                      const std::vector<std::string> &, fem::material::MaterialDatabase &)>;

using FieldRunnerMap = std::unordered_map<std::string, FieldRunner>;

FieldRunnerMap CreateFieldRunners();

int RunThermoMechanicalCoupled(const fem::frontend::ProjectConfig &config,
                               const fem::frontend::FieldConfig &thermal_field,
                               const fem::frontend::FieldConfig &mechanical_field,
                               const std::vector<std::string> &coupled_variable_names,
                               fem::material::MaterialDatabase &materials);

int RunElectroThermalCoupled(const fem::frontend::ProjectConfig &config,
                             const fem::frontend::FieldConfig &electro_field,
                             const fem::frontend::FieldConfig &thermal_field,
                             const std::vector<std::string> &coupled_variable_names,
                             fem::material::MaterialDatabase &materials);
}  // namespace fem::app
