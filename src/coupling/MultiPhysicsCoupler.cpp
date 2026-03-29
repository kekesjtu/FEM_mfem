#include "fem/coupling/MultiPhysicsCoupler.hpp"

#include "fem/coeff/Coeffmanagaer.hpp"
#include "fem/log/Logger.hpp"
#include "fem/output/ResultExporter.hpp"
#include "fem/output/SolutionTextExporter.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"

#include <cmath>

namespace fem::coupling
{

MultiPhysicsCoupler::MultiPhysicsCoupler(frontend::ProjectConfig &config) : config_(config)
{
}

void MultiPhysicsCoupler::Solve()
{
    if (config_.simulation.transient_enabled)
    {
        SolveTransient();
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
    io::ResultExporter exporter(config_.simulation.output_dir, config_.fe.GetMesh());
    if (e_solver)
        exporter.ExportScalar("electrostatic", "voltage", e_solver->GetVoltage());
    if (t_solver)
        exporter.ExportScalar("thermal", "temperature", t_solver->GetTemperature());
    if (m_solver)
        exporter.ExportVector("mechanical", "displacement", m_solver->GetDisplacement(), "u");

    logger->info("Results exported to: {}", config_.simulation.output_dir);
    logger->info("Multi-physics solve complete.");
}

// ----------------------------------------------------------------
// Transient solve — backward-Euler with Picard per step
// ----------------------------------------------------------------
void MultiPhysicsCoupler::SolveTransient()
{
    auto logger = fem::log::Get();
    const auto &simulation_config = config_.simulation;
    const bool has_e = config_.HasElectricField();
    const bool has_t = config_.HasThermalField();
    const bool has_m = config_.HasMechanicalField();

    logger->info("=== Transient solve: t=[{}, {}], dt={}, output_interval={} ===",
                 simulation_config.transient_t_start, simulation_config.transient_t_end,
                 simulation_config.transient_dt, simulation_config.transient_output_interval);

    // --- Create solvers ---
    std::unique_ptr<physics::ElectrostaticFieldSolver> e_solver;
    std::unique_ptr<physics::ThermalFieldSolver> t_solver;
    std::unique_ptr<physics::MechanicalFieldSolver> m_solver;

    if (has_e)
        e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
    if (has_t)
        t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);
    if (has_m)
        m_solver = std::make_unique<physics::MechanicalFieldSolver>(config_);

    // --- Set up coupling coefficients ---
    std::unique_ptr<coeff::ExpressionCoefficient> sigma_coeff;
    if (has_e && has_t)
    {
        sigma_coeff = std::make_unique<coeff::ExpressionCoefficient>(
            "electrical_conductivity", config_.materials, config_.fe.GetMesh(),
            &t_solver->GetTemperature());
        t_solver->SetVoltageField(&e_solver->GetVoltage());
        t_solver->SetElectricalConductivity(sigma_coeff.get());
        e_solver->SetTemperatureField(&t_solver->GetTemperature());
    }

    // T_old for backward Euler
    mfem::GridFunction T_old(&config_.fe.GetScalarFESpace());
    T_old = 293.15;

    // --- Snapshot collection ---
    std::vector<TransientSnapshot> snapshots;
    double output_interval = simulation_config.transient_output_interval;

    auto take_snapshot = [&](double t)
    {
        TransientSnapshot snap;
        snap.time = t;
        if (e_solver)
            snap.voltage = e_solver->GetVoltage();
        if (t_solver)
            snap.temperature = t_solver->GetTemperature();
        if (m_solver)
            snap.displacement = m_solver->GetDisplacement();
        snapshots.push_back(std::move(snap));
        logger->info("  Snapshot taken at t={:.4f}", t);
    };

    // --- Initial state at t=0 ---
    // Solve E with initial T (quasi-static), keep T=T_init, solve M
    if (e_solver)
        e_solver->Solve();
    if (m_solver)
    {
        if (t_solver)
            m_solver->SetTemperatureField(&t_solver->GetTemperature());
        m_solver->Solve();
    }
    take_snapshot(simulation_config.transient_t_start);
    double next_output = simulation_config.transient_t_start + output_interval;

    // --- Time stepping loop ---
    double t = simulation_config.transient_t_start;
    double dt = simulation_config.transient_dt;
    int step = 0;

    while (t < simulation_config.transient_t_end - 1e-12 * dt)
    {
        // Clamp dt to not overshoot t_end
        if (t + dt > simulation_config.transient_t_end + 1e-12)
            dt = simulation_config.transient_t_end - t;

        t += dt;
        step++;

        // Save T_old for backward Euler
        if (has_t)
            T_old = t_solver->GetTemperature();

        // Enable transient mode on thermal solver
        if (has_t)
            t_solver->EnableTransient(dt, &T_old);

        // --- Picard iteration for electro-thermal coupling ---
        int picard_iters = 0;
        bool converged = false;

        if (has_e && has_t)
        {
            mfem::Vector T_prev(t_solver->GetTemperature().Size());
            T_prev = t_solver->GetTemperature();

            for (int iter = 0; iter < simulation_config.picard_max_iterations; ++iter)
            {
                picard_iters = iter + 1;
                e_solver->Solve();
                t_solver->Solve();

                // Check Picard convergence
                mfem::Vector delta(t_solver->GetTemperature());
                delta -= T_prev;
                double diff = delta.Norml2();
                double ref = T_prev.Norml2();
                double rel_change = (ref > 0.0) ? diff / ref : diff;

                if (rel_change < simulation_config.picard_tolerance)
                {
                    converged = true;
                    break;
                }
                T_prev = t_solver->GetTemperature();
            }

            // Adaptive dt: reject step and retry with smaller dt
            if (!converged && simulation_config.adaptive_dt)
            {
                t -= dt;
                dt *= 0.5;
                if (has_t)
                    t_solver->GetTemperature() = T_old;  // restore
                logger->warn("Step {}: Picard did not converge, halving dt to {:.4e}", step, dt);
                step--;
                continue;
            }
        }
        else if (has_t)
        {
            t_solver->Solve();
            picard_iters = 1;
            converged = true;
        }
        else if (has_e)
        {
            e_solver->Solve();
            picard_iters = 1;
            converged = true;
        }

        // --- Mechanical (quasi-static with current T) ---
        if (has_m)
        {
            if (has_t)
                m_solver->SetTemperatureField(&t_solver->GetTemperature());
            m_solver->Solve();
        }

        // --- Logging ---
        logger->info("Step {}: t={:.4f}, dt={:.4e}, Picard iters={}{}", step, t, dt, picard_iters,
                     converged ? "" : " (NOT converged)");

        // --- Snapshot at output intervals ---
        if (t >= next_output - 1e-10 * output_interval)
        {
            take_snapshot(t);
            next_output += output_interval;
        }

        // --- Adaptive dt update ---
        if (simulation_config.adaptive_dt && converged)
        {
            if (picard_iters <= 3)
                dt = std::min(dt * 1.5, output_interval);
            else if (picard_iters > 6)
                dt = std::max(dt * 0.5, simulation_config.transient_dt * 0.01);
        }
    }

    logger->info("Transient solve complete. {} time steps, {} snapshots.", step, snapshots.size());

    // --- Export results ---
    ExportTransientResults(snapshots);
}

// ----------------------------------------------------------------
void MultiPhysicsCoupler::ExportTransientResults(const std::vector<TransientSnapshot> &snapshots)
{
    auto logger = fem::log::Get();
    const auto &out_dir = config_.simulation.output_dir;

    // ParaView: export last snapshot for quick visual check
    io::ResultExporter exporter(out_dir, config_.fe.GetMesh());
    if (!snapshots.empty())
    {
        const auto &last = snapshots.back();

        if (last.voltage.Size() > 0)
        {
            mfem::GridFunction gf(&config_.fe.GetScalarFESpace());
            gf = last.voltage;
            exporter.ExportScalar("electrostatic", "voltage", gf,
                                  static_cast<int>(snapshots.size()) - 1, last.time);
        }
        if (last.temperature.Size() > 0)
        {
            mfem::GridFunction gf(&config_.fe.GetScalarFESpace());
            gf = last.temperature;
            exporter.ExportScalar("thermal", "temperature", gf,
                                  static_cast<int>(snapshots.size()) - 1, last.time);
        }
        if (last.displacement.Size() > 0)
        {
            mfem::GridFunction gf(&config_.fe.GetVectorFESpace());
            gf = last.displacement;
            exporter.ExportVector("mechanical", "displacement", gf, "u",
                                  static_cast<int>(snapshots.size()) - 1, last.time);
        }
    }

    // Transient text export (serial only — parallel uses ParaView)
    if (!config_.fe.IsParallel())
    {
        const std::string txt_path = out_dir + "/transient.txt";
        std::vector<double> times;
        std::vector<mfem::Vector> voltages, temperatures, displacements;

        for (const auto &snap : snapshots)
        {
            times.push_back(snap.time);
            voltages.push_back(snap.voltage);
            temperatures.push_back(snap.temperature);
            displacements.push_back(snap.displacement);
        }

        post::SolutionTextExporter::ExportTransientNodalTxt(
            txt_path, config_.fe.GetMesh(), config_.fe.GetScalarFESpace(),
            config_.fe.GetVectorFESpace(), times, voltages, temperatures, displacements);
    }

    logger->info("Transient results exported to: {}", out_dir);
}

void MultiPhysicsCoupler::SolveElectroThermalOneWay(
    std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
    std::unique_ptr<physics::ThermalFieldSolver> &t_solver)
{
    auto logger = fem::log::Get();
    logger->info("--- Electro-thermal one-way coupling (no iteration) ---");

    // 1. Solve electrostatic field
    e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
    e_solver->Solve();

    // 2. Build sigma coefficient for Joule heating
    coeff::ExpressionCoefficient sigma_coeff("electrical_conductivity", config_.materials,
                                             config_.fe.GetMesh(), nullptr);

    // 3. Solve thermal field with Joule heating source
    t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);
    t_solver->SetVoltageField(&e_solver->GetVoltage());
    t_solver->SetElectricalConductivity(&sigma_coeff);
    t_solver->Solve();
}

void MultiPhysicsCoupler::SolveElectroThermalPicard(
    std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
    std::unique_ptr<physics::ThermalFieldSolver> &t_solver)
{
    auto logger = fem::log::Get();
    const auto &sim = config_.simulation;
    logger->info("--- Electro-thermal Picard iteration (max_iter={}, tol={:.2e}) ---",
                 sim.picard_max_iterations, sim.picard_tolerance);

    e_solver = std::make_unique<physics::ElectrostaticFieldSolver>(config_);
    t_solver = std::make_unique<physics::ThermalFieldSolver>(config_);

    // sigma depends on T — will be updated each iteration
    coeff::ExpressionCoefficient sigma_coeff("electrical_conductivity", config_.materials,
                                             config_.fe.GetMesh(), &t_solver->GetTemperature());

    t_solver->SetVoltageField(&e_solver->GetVoltage());
    t_solver->SetElectricalConductivity(&sigma_coeff);
    e_solver->SetTemperatureField(&t_solver->GetTemperature());

    mfem::Vector T_prev(t_solver->GetTemperature().Size());
    T_prev = 293.15;

    for (int iter = 0; iter < sim.picard_max_iterations; ++iter)
    {
        logger->info("Picard iteration {}/{}", iter + 1, sim.picard_max_iterations);

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

        if (rel_change < sim.picard_tolerance)
        {
            logger->info("Picard converged after {} iterations.", iter + 1);
            break;
        }

        // Relaxation
        if (std::abs(sim.picard_relaxation - 1.0) > 1e-14)
        {
            mfem::Vector &T_cur = t_solver->GetTemperature();
            for (int i = 0; i < T_cur.Size(); ++i)
            {
                T_cur(i) =
                    sim.picard_relaxation * T_cur(i) + (1.0 - sim.picard_relaxation) * T_prev(i);
            }
        }

        T_prev = t_solver->GetTemperature();
    }
}

}  // namespace fem::coupling
