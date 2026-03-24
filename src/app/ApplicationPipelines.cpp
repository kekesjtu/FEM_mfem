#include "fem/app/ApplicationPipelines.hpp"

#include "fem/app/FieldRunners.hpp"

#include "mfem.hpp"
#include "spdlog/spdlog.h"

#include <stdexcept>

namespace fem::app
{
std::vector<std::string> BuildCoupledVariableNames(const fem::frontend::ProjectConfig &config)
{
    std::vector<std::string> coupled_variable_names;
    coupled_variable_names.reserve(config.fields.size());
    for (const auto &configured_field : config.fields)
    {
        coupled_variable_names.push_back(configured_field.name);
    }
    return coupled_variable_names;
}

fem::material::MaterialDatabase BuildMaterialDatabase(
    const fem::frontend::ProjectConfig &config,
    const std::vector<std::string> &coupled_variable_names)
{
    fem::material::MaterialDatabase materials;
    for (const auto &material : config.materials)
    {
        materials.Add(material, coupled_variable_names);
    }
    return materials;
}

FieldSelection SelectFieldsByType(const fem::frontend::ProjectConfig &config)
{
    FieldSelection selection;
    for (const auto &configured_field : config.fields)
    {
        if (configured_field.type == "electrostatic" && !selection.electrostatic)
        {
            selection.electrostatic = &configured_field;
        }
        else if (configured_field.type == "thermal" && !selection.thermal)
        {
            selection.thermal = &configured_field;
        }
        else if (configured_field.type == "mechanical" && !selection.mechanical)
        {
            selection.mechanical = &configured_field;
        }
    }
    return selection;
}

int RunCoupledPipelineByPhysicsType(const fem::frontend::ProjectConfig &config,
                                    const FieldSelection &selection,
                                    const std::vector<std::string> &coupled_variable_names,
                                    fem::material::MaterialDatabase &materials)
{
#if defined(MFEM_USE_MPI)
    int rank = mfem::Mpi::WorldRank();
#else
    int rank = 0;
#endif

    if (selection.thermal && selection.mechanical)
    {
        if (rank == 0)
        {
            spdlog::info("Running thermo-mechanical coupled analysis: thermal='{}', mechanical='{}'.",
                         selection.thermal->name, selection.mechanical->name);
        }
        const int rc = RunThermoMechanicalCoupled(config, *selection.thermal, *selection.mechanical,
                                                  coupled_variable_names, materials);
        if (rank == 0)
        {
            spdlog::info("Finished coupled analysis and wrote output to '{}'.",
                         config.simulation.output_dir);
        }
        return rc;
    }

    if (selection.electrostatic && selection.thermal && !selection.mechanical)
    {
        if (rank == 0)
        {
            spdlog::info(
                "Running electro-thermal one-way coupled analysis: electrostatic='{}', thermal='{}'.",
                selection.electrostatic->name, selection.thermal->name);
        }
        const int rc =
            RunElectroThermalCoupled(config, *selection.electrostatic, *selection.thermal,
                                     coupled_variable_names, materials);
        if (rank == 0)
        {
            spdlog::info("Finished coupled analysis and wrote output to '{}'.",
                         config.simulation.output_dir);
        }
        return rc;
    }

    return -1;
}

int RunSingleFieldPipelineByPhysicsType(const fem::frontend::ProjectConfig &config,
                                        const std::vector<std::string> &coupled_variable_names,
                                        fem::material::MaterialDatabase &materials)
{
    const FieldRunnerMap runners = CreateFieldRunners();
    int final_rc = 0;

    for (const auto &configured_field : config.fields)
    {
        spdlog::info("Running field '{}' of type '{}'.", configured_field.name,
                     configured_field.type);

        const auto runner_it = runners.find(configured_field.type);
        if (runner_it == runners.end())
        {
            throw std::runtime_error("Unsupported field type: " + configured_field.type);
        }

        final_rc = runner_it->second(config, configured_field, coupled_variable_names, materials);
        if (final_rc != 0)
        {
            spdlog::error("Field '{}' failed with code {}.", configured_field.name, final_rc);
            return final_rc;
        }

        spdlog::info("Finished field '{}' and wrote ParaView output to '{}'.",
                     configured_field.name, config.simulation.output_dir);
    }

    return final_rc;
}
}  // namespace fem::app
