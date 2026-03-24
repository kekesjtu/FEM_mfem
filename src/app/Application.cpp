#include "fem/app/Application.hpp"
#include "fem/app/ApplicationPipelines.hpp"

#include "fem/frontend/ConfigLoader.hpp"
#include "fem/log/Logger.hpp"

#include "mfem.hpp"

#include "spdlog/spdlog.h"

#include <stdexcept>

namespace fem::app
{
int Application::Run(const std::string &config_path)
{
    const auto config = fem::frontend::ConfigLoader::LoadFromFile(config_path);

#if defined(MFEM_USE_MPI)
    if (config.simulation.assembly_mode == "parallel" && !mfem::Mpi::IsInitialized())
    {
        mfem::Mpi::Init();
    }
#endif

    fem::log::Init(config.simulation.log_level);

#if defined(MFEM_USE_MPI)
    int rank = mfem::Mpi::WorldRank();
#else
    int rank = 0;
#endif

    // Only rank 0 prints initialization logs to avoid duplication
    if (rank == 0)
    {
        spdlog::info("Loaded configuration: {}", config_path);
        spdlog::info("Mesh: {}, order: {}, refinement: {}", config.simulation.mesh_path,
                     config.simulation.order, config.simulation.uniform_refinement_levels);
        spdlog::info("Linear solver: {}", config.simulation.solver);
        spdlog::info("Assembly mode: {}", config.simulation.assembly_mode);

        if (config.simulation.assembly_mode == "parallel")
        {
#if defined(MFEM_USE_MPI)
            spdlog::info("MPI world size: {}, rank: {}", mfem::Mpi::WorldSize(),
                         mfem::Mpi::WorldRank());
#else
            throw std::runtime_error(
                "assembly_mode='parallel' requested, but MFEM was built without MPI support.");
#endif
        }
    }

    if (config.simulation.assembly_mode == "parallel" && rank != 0)
    {
#if defined(MFEM_USE_MPI)
        // Verify MPI is initialized on non-zero ranks
        if (!mfem::Mpi::IsInitialized())
        {
            throw std::runtime_error("MPI not initialized on non-zero rank.");
        }
#endif
    }

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
