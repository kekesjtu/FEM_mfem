#include "fem/steady/SteadySolver.hpp"

#include "fem/coupling/PicardCoupler.hpp"
#include "fem/log/Logger.hpp"
#include "fem/output/ResultExporter.hpp"
#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"

#include <mpi.h>

namespace fem::steady
{

SteadySolver::SteadySolver(frontend::ProjectConfig& config) : config_(config)
{
}

void SteadySolver::SolveSteady()
{
    auto logger = fem::log::Get();
    const bool has_e = config_.HasElectricField();
    const bool has_t = config_.HasThermalField();
    const bool has_m = config_.HasMechanicalField();
    const bool is_picard = config_.NeedsPicardIteration();

    std::unique_ptr<physics::ElectrostaticFieldSolver> e_solver;
    std::unique_ptr<physics::ThermalFieldSolver> t_solver;
    std::unique_ptr<physics::MechanicalFieldSolver> m_solver;

    // Phase 1: Electro-thermal coupling
    if (has_e && has_t && is_picard)
    {
        logger->info("--- Electro-thermal Picard iteration ---");

        e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
        t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);

        e_solver->SetTemperatureField(&t_solver->GetTemperature());
        t_solver->SetVoltageField(&e_solver->GetVoltage());
        t_solver->SetElectricalConductivity(e_solver->GetSigma());
        t_solver->CacheMatrices();

        coupling::PicardCoupler coupler(config_.simulation);
        auto result = coupler.RunPicard(e_solver.get(), t_solver.get(), has_e, has_t,
                                        coupling::PicardSolveMode::Steady);

        if (result.converged)
            logger->info("Picard converged after {} iterations.", result.iterations);
        else
            logger->warn("Picard did NOT converge after {} iterations.", result.iterations);
    }
    else if (has_e && has_t)
    {
        logger->info("--- Electro-thermal one-way coupling (no iteration) ---");

        e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
        t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);

        e_solver->SetTemperatureField(&t_solver->GetTemperature());
        e_solver->Solve();

        t_solver->SetVoltageField(&e_solver->GetVoltage());
        t_solver->SetElectricalConductivity(e_solver->GetSigma());
        t_solver->CacheMatrices();
        t_solver->SolveSteady();
    }
    else if (has_e)
    {
        e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
        e_solver->Solve();
    }
    else if (has_t)
    {
        t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);
        t_solver->CacheMatrices();
        t_solver->SolveSteady();
    }

    // Phase 2: Mechanical (optionally coupled with thermal)
    if (has_m)
    {
        m_solver = std::make_unique<physics::MechanicalFieldSolver>(config_);
        if (has_t && t_solver)
        {
            logger->debug("Thermo-mechanical coupling");
            m_solver->SetTemperatureField(&t_solver->GetTemperature());
        }
        m_solver->Solve();
    }

    // Phase 3: Export results
    MPI_Barrier(MPI_COMM_WORLD);
    output::ResultExporter exporter(config_.simulation.output_dir, config_.fe.GetMesh());
    if (e_solver)
        exporter.ExportScalar("electrostatic", "voltage", e_solver->GetVoltage());
    if (t_solver)
        exporter.ExportScalar("thermal", "temperature", t_solver->GetTemperature());
    if (m_solver)
        exporter.ExportVector("mechanical", "displacement", m_solver->GetDisplacement(), "u");

    logger->info("Results exported to: {}", config_.simulation.output_dir);
    logger->info("Steady-state solve complete.");
}

}  // namespace fem::steady
