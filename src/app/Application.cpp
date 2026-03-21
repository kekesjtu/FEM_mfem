#include "fem/app/Application.hpp"
#include "fem/app/ApplicationPipelines.hpp"

#include "fem/frontend/ConfigLoader.hpp"
#include "fem/log/Logger.hpp"

#include "spdlog/spdlog.h"

#include <stdexcept>

namespace fem::app
{
int Application::Run(const std::string &config_path)
{
    const auto config = fem::frontend::ConfigLoader::LoadFromFile(config_path);
    fem::log::Init(config.simulation.log_level);

    spdlog::info("Loaded configuration: {}", config_path);
    spdlog::info("Mesh: {}, order: {}, refinement: {}", config.simulation.mesh_path,
                 config.simulation.order, config.simulation.uniform_refinement_levels);
    spdlog::info("Linear solver: {}", config.simulation.solver);

    if (config.fields.empty())
    {
        throw std::runtime_error("No field definitions in config.");
    }

    const std::vector<std::string> coupled_variable_names = BuildCoupledVariableNames(config);
    fem::material::MaterialDatabase materials =
        BuildMaterialDatabase(config, coupled_variable_names);
    const FieldSelection selection = SelectFieldsByType(config);

    const int coupled_rc =
        RunCoupledPipelineByPhysicsType(config, selection, coupled_variable_names, materials);
    if (coupled_rc >= 0)
    {
        return coupled_rc;
    }

    return RunSingleFieldPipelineByPhysicsType(config, coupled_variable_names, materials);
}
}  // namespace fem::app
