#include "app/Application.hpp"
#include "app/PythonTool.hpp"

#include "frontend/ConfigLoader.hpp"
#include "log/Logger.hpp"
#include "steady/SteadySolver.hpp"
#include "transient/TransientSolver.hpp"

#include <mpi.h>
#include <exception>
#include <filesystem>
#include <string>

namespace fem::app
{

namespace
{
void RunComparisonTool(const std::string &script, const std::string &mine_file,
                       const std::string &reference_file, const std::string &extra_args)
{
    auto logger = fem::log::Get();
    try
    {
        const auto python = ResolvePython(std::filesystem::current_path());
        std::vector<std::string> arguments{script, "--mine", mine_file, "--comsol", reference_file};
        const auto extra = SplitToolArguments(extra_args);
        arguments.insert(arguments.end(), extra.begin(), extra.end());
        logger->info("Python interpreter: {}", python.executable.u8string());
        const int status = RunPython(python, arguments);
        if (status != 0)
            logger->warn("Comparison script returned non-zero exit code: {}", status);
    }
    catch (const std::exception &error)
    {
        logger->warn("Comparison could not run: {}", error.what());
    }
}

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

    logger->info("=== Running COMSOL comparison ===");
    RunComparisonTool("tools/compare_solution_txt.py", mine_file,
                      sim.comsol_reference_path, sim.compare_args);
}

void RunComsolCurveComparison(const frontend::SimulationConfig &sim)
{
    if (sim.curve_reference_path.empty())
        return;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0)
        return;

    auto logger = fem::log::Get();

    if (!sim.transient_enabled)
    {
        logger->warn("Curve comparison skipped: transient_enabled is false.");
        return;
    }

    const std::string mine_file = sim.output_dir + "/transient.txt";
    if (!std::filesystem::exists(mine_file))
    {
        logger->warn("Curve comparison skipped: output file not found: {}", mine_file);
        return;
    }
    if (!std::filesystem::exists(sim.curve_reference_path))
    {
        logger->warn("Curve comparison skipped: COMSOL reference not found: {}",
                     sim.curve_reference_path);
        return;
    }

    logger->info("=== Running COMSOL time-curve comparison ===");
    RunComparisonTool("tools/compare_time_curves.py", mine_file,
                      sim.curve_reference_path, sim.curve_compare_args);
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
    RunComsolCurveComparison(config.simulation);

    logger->info("All done.");
    return 0;
}

}  // namespace fem::app
