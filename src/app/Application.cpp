#include "fem/app/Application.hpp"

#include "fem/frontend/ConfigLoader.hpp"
#include "fem/log/Logger.hpp"
#include "fem/steady/SteadySolver.hpp"
#include "fem/transient/TransientSolver.hpp"

#include <mpi.h>
#include <cstdlib>
#include <filesystem>
#include <string>

namespace fem::app
{

namespace
{
void RunComsolComparison(const frontend::SimulationConfig &sim)
{
    if (sim.comsol_reference_path.empty())
        return;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0)
        return;

    auto logger = fem::log::Get();

    // Determine the mine file path
    std::string mine_file;
    if (sim.compare_args.find("--transient") != std::string::npos)
    {
        mine_file = sim.output_dir + "/transient.txt";
    }
    else
    {
        // Guess from COMSOL reference filename
        const auto ref = sim.comsol_reference_path;
        if (ref.find("displacement") != std::string::npos)
            mine_file = sim.output_dir + "/displacement.txt";
        else
            mine_file = sim.output_dir + "/temperature.txt";
    }

    if (!std::filesystem::exists(mine_file))
    {
        logger->warn("Compare skipped: output file not found: {}", mine_file);
        return;
    }
    if (!std::filesystem::exists(sim.comsol_reference_path))
    {
        logger->warn("Compare skipped: COMSOL reference not found: {}", sim.comsol_reference_path);
        return;
    }

    std::string cmd = "python3 tools/compare_solution_txt.py --mine \"" + mine_file +
                      "\" --comsol \"" + sim.comsol_reference_path + "\"";
    if (!sim.compare_args.empty())
        cmd += " " + sim.compare_args;

    logger->info("=== Running COMSOL comparison ===");
    logger->info("Command: {}", cmd);
    int ret = std::system(cmd.c_str());
    if (ret != 0)
        logger->warn("Comparison script returned non-zero exit code: {}", ret);
}
}  // namespace

int Application::Run(const std::string &config_path)
{
    auto config = frontend::ConfigLoader::LoadFromFile(config_path);
    auto logger = fem::log::Get();

    logger->info("====================================");
    logger->info("  FEM Multi-Physics Solver");
    logger->info("====================================");

    if (config.simulation.transient_enabled)
    {
        transient::TransientSolver transient_solver(config);
        transient_solver.SolveTransient();
    }
    else
    {
        steady::SteadySolver steady_solver(config);
        steady_solver.SolveSteady();
    }

    RunComsolComparison(config.simulation);

    logger->info("All done.");
    return 0;
}

}  // namespace fem::app
