#include "fem/coupling/MultiPhysicsCoupler.hpp"

#include "fem/coeff/Coeffmanagaer.hpp"
#include "fem/log/Logger.hpp"
#include "fem/output/ResultExporter.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/transient/TransientSolver.hpp"

namespace fem::coupling
{

MultiPhysicsCoupler::MultiPhysicsCoupler(frontend::ProjectConfig &config) : config_(config)
{
}

void MultiPhysicsCoupler::Solve()
{
    if (config_.simulation.transient_enabled)
    {
        transient::TransientSolver ts(config_);
        ts.Run();
    }
    else
    {
        SolveSteadyState();
    }
}

// ----------------------------------------------------------------
// Steady-state solve (original path, unchanged semantics)
// ----------------------------------------------------------------
void MultiPhysicsCoupler::SolveSteadyState()
{
    auto logger = fem::log::Get();
    const bool has_e = config_.HasElectricField();
    const bool has_t = config_.HasThermalField();
    const bool has_m = config_.HasMechanicalField();
    const bool picard = config_.NeedsPicardIteration();

    std::unique_ptr<physics::ElectrostaticFieldSolver> e_solver;
    std::unique_ptr<physics::ThermalFieldSolver> t_solver;
    std::unique_ptr<physics::MechanicalFieldSolver> m_solver;

    // Phase 1: Electro-thermal coupling
    if (has_e && has_t)
    {
        if (picard)
            SolveElectroThermalPicard(e_solver, t_solver);
        else
            SolveElectroThermalOneWay(e_solver, t_solver);
    }
    else if (has_e)
    {
        e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
        e_solver->Solve();
    }
    else if (has_t)
    {
        t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);
        t_solver->Solve();
    }

    // Phase 2: Mechanical (optionally coupled with thermal)
    if (has_m)
    {
        m_solver = std::make_unique<physics::MechanicalFieldSolver>(config_);
        if (has_t && t_solver)
        {
            logger->info("--- Thermo-mechanical coupling ---");
            m_solver->SetTemperatureField(&t_solver->GetTemperature());
        }
        m_solver->Solve();
    }

    // Phase 3: Export all results via unified exporter
    output::ResultExporter exporter(config_.simulation.output_dir, config_.fe.GetMesh());
    if (e_solver)
        exporter.ExportScalar("electrostatic", "voltage", e_solver->GetVoltage());
    if (t_solver)
        exporter.ExportScalar("thermal", "temperature", t_solver->GetTemperature());
    if (m_solver)
        exporter.ExportVector("mechanical", "displacement", m_solver->GetDisplacement(), "u");

    logger->info("Results exported to: {}", config_.simulation.output_dir);
    logger->info("Multi-physics solve complete.");
}

void MultiPhysicsCoupler::SolveElectroThermalOneWay(
    std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
    std::unique_ptr<physics::ThermalFieldSolver> &t_solver)
{
    auto logger = fem::log::Get();
    logger->info("--- Electro-thermal one-way coupling (no iteration) ---");
    
    e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
    t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);

    // 2. Build sigma coefficient for Joule heating
    coeff::ExpressionCoefficient sigma_coeff("electrical_conductivity", config_.materials,
                                             config_.fe.GetMesh(), nullptr);
    // 1. Solve electrostatic field
    e_solver->SetTemperatureField(&t_solver->GetTemperature());
    e_solver->SetElectricalConductivity(&sigma_coeff);
    e_solver->Solve();

    // 3. Solve thermal field with Joule heating source
    t_solver->SetVoltageField(&e_solver->GetVoltage());
    t_solver->SetElectricalConductivity(&sigma_coeff);
    t_solver->Solve();
}

void MultiPhysicsCoupler::SolveElectroThermalPicard(
    std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
    std::unique_ptr<physics::ThermalFieldSolver> &t_solver)
{
    auto logger = fem::log::Get();
    const auto &simulation_config = config_.simulation;
    const double T_init = config_.thermal_field.initial_temperature;
    logger->info("--- Electro-thermal Picard iteration (max_iter={}, tol={:.2e}) ---",
                 simulation_config.picard_max_iterations, simulation_config.picard_tolerance);

    e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
    t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);

    // sigma depends on T — will be updated each iteration
    coeff::ExpressionCoefficient sigma_coeff("electrical_conductivity", config_.materials,
                                             config_.fe.GetMesh(), &t_solver->GetTemperature());

    e_solver->SetTemperatureField(&t_solver->GetTemperature());
    // Share the same sigma coefficient with e_solver so the expression
    // "5.998e7/(1+0.00393*(T-293.15))" is evaluated only once per quadrature point.
    e_solver->SetElectricalConductivity(&sigma_coeff);

    t_solver->SetVoltageField(&e_solver->GetVoltage());
    t_solver->SetElectricalConductivity(&sigma_coeff);

    mfem::Vector T_prev(t_solver->GetTemperature().Size());
    T_prev = T_init;

    for (int iter = 0; iter < simulation_config.picard_max_iterations; ++iter)
    {
        logger->info("Picard iteration {}/{}", iter + 1, simulation_config.picard_max_iterations);

        // 1. Solve electric field with current T
        e_solver->Solve();

        // 2. Solve thermal field with Joule heating from V
        t_solver->Solve();

        // 3. Check convergence
        mfem::Vector delta_T(t_solver->GetTemperature());
        delta_T -= T_prev;
        double diff_norm = delta_T.Norml2();
        double ref_norm = T_prev.Norml2();
        double rel_change = (ref_norm > 0.0) ? diff_norm / ref_norm : diff_norm;

        logger->info("  Picard rel_change = {:.6e}", rel_change);

        if (rel_change < simulation_config.picard_tolerance)
        {
            logger->info("Picard converged after {} iterations.", iter + 1);
            break;
        }

        // Relaxation
        if (std::abs(simulation_config.picard_relaxation - 1.0) > 1e-14)
        {
            mfem::Vector &T_cur = t_solver->GetTemperature();
            for (int i = 0; i < T_cur.Size(); ++i)
            {
                T_cur(i) = simulation_config.picard_relaxation * T_cur(i) +
                           (1.0 - simulation_config.picard_relaxation) * T_prev(i);
            }
        }

        T_prev = t_solver->GetTemperature();
    }
}

}  // namespace fem::coupling
