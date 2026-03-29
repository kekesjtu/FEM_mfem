#include "fem/app/Application.hpp"

#include "fem/coupling/MultiPhysicsCoupler.hpp"
#include "fem/frontend/ConfigLoader.hpp"
#include "fem/log/Logger.hpp"

namespace fem::app
{

int Application::Run(const std::string &config_path)
{
    auto config = frontend::ConfigLoader::LoadFromFile(config_path);
    auto logger = fem::log::Get();

    logger->info("====================================");
    logger->info("  FEM Multi-Physics Solver");
    logger->info("====================================");
    logger->info("Config: {}", config_path);
    logger->info("Mesh: {}", config.simulation.mesh_path);
    logger->info("Order: {}", config.simulation.order);
    logger->info("Solver: {}", config.simulation.solver);

    coupling::MultiPhysicsCoupler coupler(config);
    coupler.Solve();

    logger->info("All done.");
    return 0;
}

}  // namespace fem::app
