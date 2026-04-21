#include "fem/transient/TransientSolver.hpp"

#include "fem/coeff/CoefficientManager.hpp"
#include "fem/coupling/PicardCoupler.hpp"
#include "fem/log/Logger.hpp"
#include "fem/output/ResultExporter.hpp"
#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"
#include "fem/transient/TransientUtils.hpp"

#include <algorithm>
#include <cmath>

namespace fem::transient
{

TransientSolver::TransientSolver(frontend::ProjectConfig &config) : config_(config)
{
}

// ----------------------------------------------------------------
// Transient solve — predictor-corrector adaptive dt with order selection
// Phase 1: Electro-thermal coupling (all time steps)
// Phase 2: Mechanical batch solve (cached stiffness matrix, RHS-only)
// ----------------------------------------------------------------
void TransientSolver::SolveTransient()
{
    auto logger = fem::log::Get();
    const auto &sim = config_.simulation;
    const bool has_e = config_.HasElectricField();
    const bool has_t = config_.HasThermalField();
    const bool has_m = config_.HasMechanicalField();
    const double dt_max = (sim.dt_max > 0.0) ? sim.dt_max : sim.transient_output_interval;
    const double dt_min = sim.dt_min;
    const double T_init = config_.thermal_field.initial_temperature;

    logger->info("=== Transient solve: t=[{}, {}], dt={}, interval={} ===", sim.transient_t_start,
                 sim.transient_t_end, sim.transient_dt, sim.transient_output_interval);
    if (sim.adaptive_dt)
    {
        logger->info("  Adaptive PC: reltol={:.2e} abstol={:.2e} dt=[{:.2e},{:.2e}]",
                     sim.adaptive_reltol, sim.adaptive_abstol, dt_min, dt_max);
    }

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

    // --- Coupling coefficients ---
    std::unique_ptr<coeff::ExpressionCoefficient> sigma_coeff;
    if (has_e && has_t)
    {
        sigma_coeff = std::make_unique<coeff::ExpressionCoefficient>(
            "electrical_conductivity", config_.materials, config_.fe.GetMesh(),
            &t_solver->GetTemperature());
        e_solver->SetTemperatureField(&t_solver->GetTemperature());
        e_solver->SetElectricalConductivity(sigma_coeff.get());
        t_solver->SetVoltageField(&e_solver->GetVoltage());
        t_solver->SetElectricalConductivity(sigma_coeff.get());
    }

    // T_old for backward Euler
    mfem::ParGridFunction T_old(&config_.fe.GetScalarFESpace());
    T_old = T_init;

    // History for predictor-corrector error estimation and adaptive BDF order
    const int N_true = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector T_nm1_true(N_true);  // T^{n-1}
    mfem::Vector T_nm2_true(N_true);  // T^{n-2}
    T_nm1_true = 0.0;
    T_nm2_true = 0.0;
    double dt_prev = 0.0;
    double dt_prev2 = 0.0;
    int history_depth = 0;  // 0=cold start, 1=linear pred, 2+=quad pred+BDF2

    // ================================================================
    // Phase 1: Electro-thermal time stepping
    // ================================================================
    logger->info("Phase 1: Electro-thermal time stepping");

    const double output_interval = sim.transient_output_interval;
    double next_output_time = sim.transient_t_start + output_interval;
    std::vector<TransientSnapshot> output_snapshots;

    // Cache matrices + initial solve
    if (has_t)
        t_solver->CacheMatrices();
    if (e_solver)
    {
        e_solver->UpdateBoundaryConditions(sim.transient_t_start);
        e_solver->Solve();
    }

    // Initial snapshot
    FieldSnapshot prev_snap =
        SaveFieldSnapshot(sim.transient_t_start, e_solver.get(), t_solver.get());
    {
        TransientSnapshot snap;
        snap.time = sim.transient_t_start;
        snap.voltage = prev_snap.voltage;
        snap.temperature = prev_snap.temperature;
        output_snapshots.push_back(std::move(snap));
    }

    // --- Initial dt (consistent initialization) ---
    double t = sim.transient_t_start;
    double dt = sim.adaptive_dt ? dt_min : sim.transient_dt;

    if (sim.adaptive_dt && has_t)
    {
        double rate = t_solver->ComputeInitialRate();
        // Hairer-Wanner style: dt ≈ 0.01 * d0/d1, with d1 = max(rate, d0/dt_max)
        // When rate ≈ 0 (time-varying BC): dt_init ≈ 0.01*dt_max (conservative)
        // When rate is large: dt_init ≈ 0.01*d0/rate (small, safe)
        double d0 = sim.adaptive_abstol + sim.adaptive_reltol * std::abs(T_init);
        double d1 = std::max(rate, d0 / dt_max);
        double dt_init = std::clamp(0.01 * d0 / d1, dt_min, dt_max);
        dt = dt_init;
        logger->info("Consistent init: max |dT/dt|_0={:.4e} K/s, dt_init={:.4e}", rate, dt_init);
    }

    int step = 0;
    int rejected_steps = 0;

    // Picard coupler for E-T coupling at each time step
    coupling::PicardCoupler picard_coupler(sim);

    // --- Time stepping loop ---
    // Run until all output times are covered (may overshoot t_end)
    const double last_output_time =
        sim.transient_t_start +
        std::floor((sim.transient_t_end - sim.transient_t_start) / output_interval + 0.5) *
            output_interval;
    while (next_output_time <= last_output_time + 1e-12 * output_interval)
    {
        dt = std::clamp(dt, dt_min, dt_max);

        if (has_t)
            T_old = t_solver->GetTemperature();
        if (has_t)
            t_solver->EnableTransient(dt, &T_old);
        if (e_solver)
            e_solver->UpdateBoundaryConditions(t + dt);

        // --- Picard coupling ---
        auto picard = picard_coupler.RunPicard(e_solver.get(), t_solver.get(), has_e, has_t,
                                               coupling::PicardSolveMode::BackwardEuler);

        if (!picard.converged && sim.adaptive_dt)
        {
            if (dt <= dt_min * (1.0 + 1e-10))
            {
                logger->warn("Step force-accepted at dt_min={:.4e}: Picard not converged", dt);
            }
            else
            {
                if (has_t)
                    t_solver->GetTemperature() = T_old;
                dt *= 0.5;
                rejected_steps++;
                logger->warn("Step rejected: Picard not converged, dt->{:.4e}", dt);
                continue;
            }
        }

        // --- Predictor-corrector error estimation with adaptive BDF order ---
        double wrms_err = 0.0;
        double dt_new = dt;
        int accepted_order = 1;
        mfem::Vector T_n_true;  // T^n in true DOFs (shared by error est & history)

        if (has_t && sim.adaptive_dt)
        {
            if (history_depth == 0)
            {
                // Cold start: no predictor available, no error estimation
                dt_new = std::min(dt * sim.eta_max, dt_max);
                accepted_order = 1;
            }
            else
            {
                T_n_true.SetSize(N_true);
                T_old.GetTrueDofs(T_n_true);
                mfem::Vector T_be_true(N_true);
                t_solver->GetTemperature().GetTrueDofs(T_be_true);

                // --- Order 1: linear predictor + BE corrector ---
                mfem::Vector T_pred1(N_true);
                {
                    double r = dt / dt_prev;
                    add(1.0 + r, T_n_true, -r, T_nm1_true, T_pred1);
                }
                auto err1 = ComputeAdaptiveError(T_pred1, T_be_true, dt, sim.adaptive_reltol,
                                                 sim.adaptive_abstol, sim.eta_safety, sim.eta_max,
                                                 sim.eta_min, dt_min, 1);
                double wrms_1 = err1.wrms;
                double dt_new_1 = err1.dt_new;

                // --- Order 2: quadratic predictor + BDF2 corrector ---
                double wrms_2 = 1e30;
                double dt_new_2 = 0.0;
                if (history_depth >= 2)
                {
                    t_solver->SolveBDF2(T_nm1_true, dt_prev);

                    mfem::Vector T_pred2_true(N_true);
                    {
                        double h1 = dt_prev, h2 = dt_prev2;
                        double L0 = dt * (dt + h1) / (h2 * (h1 + h2));
                        double L1 = -dt * (dt + h1 + h2) / (h1 * h2);
                        double L2 = (dt + h1) * (dt + h1 + h2) / (h1 * (h1 + h2));
                        add(L0, T_nm2_true, L1, T_nm1_true, T_pred2_true);
                        T_pred2_true.Add(L2, T_n_true);
                    }
                    auto err2 = ComputeAdaptiveError(
                        T_pred2_true, t_solver->GetBDF2Solution(), dt, sim.adaptive_reltol,
                        sim.adaptive_abstol, sim.eta_safety, sim.eta_max, sim.eta_min, dt_min, 2);
                    wrms_2 = err2.wrms;
                    dt_new_2 = err2.dt_new;
                }

                // --- Adaptive order selection: pick order giving largest acceptable dt ---
                bool accept_order1 = (wrms_1 <= 1.0);
                bool accept_order2 = (history_depth >= 2) && (wrms_2 <= 1.0);

                if (accept_order1 && accept_order2)
                {
                    accepted_order = (dt_new_2 >= dt_new_1) ? 2 : 1;
                    wrms_err = (accepted_order == 2) ? wrms_2 : wrms_1;
                    dt_new = (accepted_order == 2) ? dt_new_2 : dt_new_1;
                }
                else if (accept_order2)
                {
                    accepted_order = 2;
                    wrms_err = wrms_2;
                    dt_new = dt_new_2;
                }
                else if (accept_order1)
                {
                    accepted_order = 1;
                    wrms_err = wrms_1;
                    dt_new = dt_new_1;
                }
                else
                {
                    // Both rejected
                    wrms_err = wrms_1;
                    dt_new = (history_depth >= 2) ? std::max(dt_new_1, dt_new_2) : dt_new_1;

                    if (dt <= dt_min * (1.0 + 1e-10))
                    {
                        logger->warn(
                            "Step force-accepted at dt_min={:.4e}: wrms1={:.4e} wrms2={:.4e}", dt,
                            wrms_1, wrms_2);
                        accepted_order = 1;
                    }
                    else
                    {
                        t_solver->GetTemperature() = T_old;
                        rejected_steps++;
                        dt = dt_new;
                        logger->info("Step rejected: wrms1={:.4e} wrms2={:.4e}, dt->{:.4e}", wrms_1,
                                     wrms_2, dt);
                        continue;
                    }
                }

                // Apply the chosen order's solution
                if (accepted_order == 2)
                    t_solver->GetTemperature().Distribute(t_solver->GetBDF2Solution());
            }
        }

        // --- Step accepted: save history ---
        if (has_t && sim.adaptive_dt)
        {
            if (T_n_true.Size() == 0)
            {
                T_n_true.SetSize(N_true);
                T_old.GetTrueDofs(T_n_true);
            }
            T_nm2_true = T_nm1_true;
            T_nm1_true = T_n_true;
            dt_prev2 = dt_prev;
            dt_prev = dt;
            history_depth = std::min(history_depth + 1, 2);
        }

        t += dt;
        step++;

        FieldSnapshot curr_snap = SaveFieldSnapshot(t, e_solver.get(), t_solver.get());

        // Emit output snapshots at interval boundaries
        while (next_output_time <= t + 1e-12 * output_interval &&
               next_output_time <= last_output_time + 1e-12 * output_interval)
        {
            EmitOutputSnapshot(next_output_time, prev_snap, curr_snap, output_snapshots);
            next_output_time += output_interval;
        }
        prev_snap = std::move(curr_snap);

        // Step log
        if (sim.adaptive_dt && has_t)
            logger->info(
                "Step {:>3} | t={:.4f} dt={:.4e} picard={} order={} wrms={:.4e} dt_next={:.4e}",
                step, t, dt, picard.iterations, accepted_order, wrms_err, dt_new);
        else
            logger->info("Step {:>3} | t={:.4f} dt={:.4e} picard={}{}", step, t, dt,
                         picard.iterations, picard.converged ? "" : " [!]");

        if (sim.adaptive_dt)
            dt = dt_new;
    }

    logger->info("Phase 1 done: {} accepted, {} rejected, {} snapshots", step, rejected_steps,
                 output_snapshots.size());

    // ================================================================
    // Phase 2: Mechanical batch solve
    // ================================================================
    if (has_m)
    {
        logger->info("Phase 2: Mechanical batch solve ({} snapshots)", output_snapshots.size());

        if (has_t)
            m_solver->SetTemperatureField(&t_solver->GetTemperature());
        m_solver->CacheStiffnessMatrix();

        for (auto &snap : output_snapshots)
        {
            if (has_t)
                t_solver->GetTemperature() = snap.temperature;
            m_solver->SolveRHSOnly();
            snap.displacement = m_solver->GetDisplacement();
        }
        logger->info("Phase 2 done: {} mechanical solves", output_snapshots.size());
    }

    logger->info("Transient complete: {} steps, {} rejected, {} snapshots", step, rejected_steps,
                 output_snapshots.size());
    ExportTransientResults(std::move(output_snapshots));
}

// ----------------------------------------------------------------
void TransientSolver::ExportTransientResults(std::vector<TransientSnapshot> snapshots)
{
    auto logger = fem::log::Get();
    const auto &out_dir = config_.simulation.output_dir;
    output::ResultExporter exporter(out_dir, config_.fe.GetMesh());

    std::vector<double> times;
    std::vector<mfem::Vector> voltages;
    std::vector<mfem::Vector> temperatures;
    std::vector<mfem::Vector> displacements;
    times.reserve(snapshots.size());
    voltages.reserve(snapshots.size());
    temperatures.reserve(snapshots.size());
    displacements.reserve(snapshots.size());

    for (auto &snap : snapshots)
    {
        times.push_back(snap.time);
        voltages.push_back(std::move(snap.voltage));
        temperatures.push_back(std::move(snap.temperature));
        displacements.push_back(std::move(snap.displacement));
    }

    if (!snapshots.empty())
    {
        exporter.ExportTransient(config_.fe.GetScalarFESpace(), config_.fe.GetVectorFESpace(),
                                 times, voltages, temperatures, displacements);
    }

    logger->info("Transient results exported to: {}", out_dir);
}

}  // namespace fem::transient
