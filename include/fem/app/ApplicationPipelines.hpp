#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/material/MaterialDatabase.hpp"

#include <string>
#include <vector>

namespace fem::app
{
struct FieldSelection
{
    const fem::frontend::FieldConfig *electrostatic = nullptr;
    const fem::frontend::FieldConfig *thermal = nullptr;
    const fem::frontend::FieldConfig *mechanical = nullptr;
};

std::vector<std::string> BuildCoupledVariableNames(const fem::frontend::ProjectConfig &config);

fem::material::MaterialDatabase BuildMaterialDatabase(
    const fem::frontend::ProjectConfig &config,
    const std::vector<std::string> &coupled_variable_names);

FieldSelection SelectFieldsByType(const fem::frontend::ProjectConfig &config);

int RunCoupledPipelineByPhysicsType(const fem::frontend::ProjectConfig &config,
                                    const FieldSelection &selection,
                                    const std::vector<std::string> &coupled_variable_names,
                                    fem::material::MaterialDatabase &materials);

int RunSingleFieldPipelineByPhysicsType(const fem::frontend::ProjectConfig &config,
                                        const std::vector<std::string> &coupled_variable_names,
                                        fem::material::MaterialDatabase &materials);
}  // namespace fem::app
