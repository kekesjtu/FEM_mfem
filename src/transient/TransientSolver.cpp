#include "fem/transient/TransientSolver.hpp"

#include "fem/coeff/Coeffmanagaer.hpp"
#include "fem/log/Logger.hpp"
#include "fem/output/ResultExporter.hpp"
#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"

#include <cmath>
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace fem::transient
{

TransientSolver::TransientSolver(frontend::ProjectConfig &config) : config_(config)
{
}

// ----------------------------------------------------------------
// Transient solve — backward-Euler with predictor-corrector adaptive dt
// ----------------------------------------------------------------
void TransientSolver::Run()
{
    auto logger = fem::log::Get();
    const auto &sim = config_.simulation;
    const bool has_e = config_.HasElectricField();
    const bool has_t = config_.HasThermalField();
    const bool has_m = config_.HasMechanicalField();
    const double dt_max = (sim.dt_max > 0.0) ? sim.dt_max : sim.transient_output_interval;
    const double dt_min = sim.dt_min;
    const double T_init = config_.thermal_field.initial_temperature;
    // backward-Euler is order k=1
    constexpr int k = 1;

    logger->info(
        "=== Transient solve: t=[{}, {}], dt={}, output_interval={} ===", sim.transient_t_start,
        sim.transient_t_end, sim.transient_dt, sim.transient_output_interval);
    if (sim.adaptive_dt)
    {
        logger->info("  Adaptive dt: reltol={:.2e}, abstol={:.2e}, dt_range=[{:.2e}, {:.2e}]",
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

    // --- Set up coupling coefficients ---
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

    // T_old for backward Euler mass matrix
    mfem::GridFunction T_old(&config_.fe.GetScalarFESpace());
    T_old = T_init;

    // History for prediction: store T at t_{n} and t_{n-1}
    mfem::Vector T_n, T_nm1;  // T at current and previous accepted steps
    double dt_prev = 0.0;     // dt of previous accepted step

    // --- Raw snapshot collection (at actual computed times) ---
    std::vector<TransientSnapshot> raw_snapshots;

    auto take_raw_snapshot = [&](double t)
    {
        TransientSnapshot snap;
        snap.time = t;
        if (e_solver)
            snap.voltage = e_solver->GetVoltage();
        if (t_solver)
            snap.temperature = t_solver->GetTemperature();
        if (m_solver)
            snap.displacement = m_solver->GetDisplacement();
        raw_snapshots.push_back(std::move(snap));
    };

    // --- Initial state at t=0 ---
    if (e_solver)
    {
        e_solver->UpdateBoundaryConditions(sim.transient_t_start);
        e_solver->Solve();
    }
    if (m_solver)
    {
        if (t_solver)
            m_solver->SetTemperatureField(&t_solver->GetTemperature());
        m_solver->Solve();
    }
    take_raw_snapshot(sim.transient_t_start);
    logger->info("  Initial snapshot at t={:.4f}", sim.transient_t_start);

    // Initialize history
    if (has_t)
    {
        T_n.SetSize(t_solver->GetTemperature().Size());
        T_n = t_solver->GetTemperature();
        T_nm1.SetSize(T_n.Size());
        T_nm1 = T_n;  // no previous step available yet
    }

    // --- Time stepping loop ---
    double t = sim.transient_t_start;
    double dt = sim.transient_dt;
    int step = 0;
    int rejected_steps = 0;

    while (t < sim.transient_t_end - 1e-12 * dt)
    {
        // Clamp dt
        dt = std::min(dt, dt_max);
        dt = std::max(dt, dt_min);
        if (t + dt > sim.transient_t_end + 1e-12)
            dt = sim.transient_t_end - t;

        // --- Step 1: Prediction ---
        mfem::Vector T_predicted;
        bool has_prediction = false;
        if (has_t && sim.adaptive_dt)
        {
            T_predicted.SetSize(T_n.Size());
            if (dt_prev > 0.0)
            {
                // Linear extrapolation: T_pred = T_n + dt * (T_n - T_{n-1}) / dt_prev
                for (int i = 0; i < T_n.Size(); i++)
                {
                    T_predicted(i) = T_n(i) + dt * (T_n(i) - T_nm1(i)) / dt_prev;
                }
            }
            else
            {
                // Zero-order extrapolation for first step: T_pred = T_n
                T_predicted = T_n;
            }
            has_prediction = true;
        }

        // Save state for backward Euler
        if (has_t)
            T_old = t_solver->GetTemperature();

        // Enable transient mode on thermal solver
        if (has_t)
            t_solver->EnableTransient(dt, &T_old);

        // --- Step 2: Correction (Picard iteration) ---
        int picard_iters = 0;
        bool picard_converged = false;

        // Update time-varying Dirichlet BCs for the target time
        if (e_solver)
            e_solver->UpdateBoundaryConditions(t + dt);

        if (has_e && has_t)
        {
            mfem::Vector T_prev(t_solver->GetTemperature().Size());
            T_prev = t_solver->GetTemperature();

            for (int iter = 0; iter < sim.picard_max_iterations; ++iter)
            {
                picard_iters = iter + 1;
                e_solver->Solve();
                t_solver->Solve();

                mfem::Vector delta(t_solver->GetTemperature());
                delta -= T_prev;
                double diff_sq = delta * delta;
                double ref_sq = T_prev * T_prev;
#ifdef MFEM_USE_MPI
                if (config_.fe.IsParallel())
                {
                    double buf[2] = {diff_sq, ref_sq};
                    double gbuf[2];
                    MPI_Allreduce(buf, gbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    diff_sq = gbuf[0];
                    ref_sq = gbuf[1];
                }
#endif
                double diff = std::sqrt(diff_sq);
                double ref = std::sqrt(ref_sq);
                double rel_change = (ref > 0.0) ? diff / ref : diff;

                if (rel_change < sim.picard_tolerance)
                {
                    picard_converged = true;
                    break;
                }
                T_prev = t_solver->GetTemperature();
            }

            // If Picard failed, reject step
            if (!picard_converged)
            {
                if (sim.adaptive_dt)
                {
                    // Restore state and halve dt
                    if (has_t)
                        t_solver->GetTemperature() = T_old;
                    dt *= 0.5;
                    rejected_steps++;
                    logger->warn("Step rejected: Picard did not converge, halving dt to {:.4e}",
                                 dt);
                    continue;
                }
                // Non-adaptive: continue with warning
            }
        }
        else if (has_t)
        {
            t_solver->Solve();
            picard_iters = 1;
            picard_converged = true;
        }
        else if (has_e)
        {
            e_solver->Solve();
            picard_iters = 1;
            picard_converged = true;
        }

        // --- Step 3 & 4: WRMS error estimation and dt control ---
        double wrms_err = 0.0;
        double dt_new = dt;

        if (has_t && sim.adaptive_dt && has_prediction)
        {
            // Compute WRMS norm: err = sqrt(1/N * sum((e_i / (abstol + reltol*|T_i|))^2))
            const mfem::Vector &T_corrected = t_solver->GetTemperature();
            int local_N = T_corrected.Size();
            double local_sum_sq = 0.0;
            for (int i = 0; i < local_N; i++)
            {
                double e_i = std::abs(T_corrected(i) - T_predicted(i));
                double weight =
                    sim.adaptive_abstol + sim.adaptive_reltol * std::abs(T_corrected(i));
                double ratio = e_i / weight;
                local_sum_sq += ratio * ratio;
            }
            double global_sum_sq = local_sum_sq;
            long long global_N = local_N;
#ifdef MFEM_USE_MPI
            if (config_.fe.IsParallel())
            {
                MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);
                long long ln = local_N;
                MPI_Allreduce(&ln, &global_N, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            }
#endif
            wrms_err = std::sqrt(global_sum_sq / global_N);

            // Step 5: Step size control decision
            if (wrms_err <= 1.0)
            {
                // Accept step
                // Optimal dt: dt_new = dt * (1/err)^(1/(k+1)) with safety factor
                double factor =
                    (wrms_err > 1e-10) ? 0.9 * std::pow(1.0 / wrms_err, 1.0 / (k + 1)) : 2.0;
                factor = std::min(factor, 2.0);  // max growth factor
                factor = std::max(factor, 0.5);  // min factor
                dt_new = dt * factor;
            }
            else
            {
                // Reject step
                double factor = 0.9 * std::pow(1.0 / wrms_err, 1.0 / (k + 1));
                factor = std::max(factor, 0.2);  // don't shrink more than 5x
                dt_new = dt * factor;
                dt_new = std::max(dt_new, dt_min);

                // Restore state
                if (has_t)
                    t_solver->GetTemperature() = T_old;
                rejected_steps++;
                dt = dt_new;
                logger->info("Step rejected: WRMS err={:.4e} > 1, new dt={:.4e}", wrms_err, dt);
                continue;
            }
        }

        // --- Step accepted: advance time ---
        t += dt;
        step++;

        // Solve mechanical (quasi-static with current T)
        if (has_m)
        {
            if (has_t)
                m_solver->SetTemperatureField(&t_solver->GetTemperature());
            m_solver->Solve();
        }

        // --- Update history ---
        if (has_t)
        {
            T_nm1 = T_n;
            T_n = t_solver->GetTemperature();
            dt_prev = dt;
        }

        // --- Logging ---
        if (sim.adaptive_dt && has_prediction)
        {
            logger->info("Step {}: t={:.4f}, dt={:.4e}, Picard={}, WRMS={:.4e}, dt_next={:.4e}",
                         step, t, dt, picard_iters, wrms_err, dt_new);
        }
        else
        {
            logger->info("Step {}: t={:.4f}, dt={:.4e}, Picard={}{}", step, t, dt, picard_iters,
                         picard_converged ? "" : " (NOT converged)");
        }

        // --- Take raw snapshot at every accepted step ---
        take_raw_snapshot(t);

        // --- Update dt for next step ---
        if (sim.adaptive_dt)
        {
            dt = dt_new;
        }
    }

    logger->info("Transient solve complete. {} accepted steps, {} rejected, {} raw snapshots.",
                 step, rejected_steps, raw_snapshots.size());

    // --- Interpolate to uniform output grid and export ---
    if (sim.adaptive_dt)
    {
        auto uniform_snapshots =
            InterpolateSnapshots(raw_snapshots, sim.transient_t_start, sim.transient_t_end,
                                 sim.transient_output_interval);
        logger->info("Interpolated to {} uniform snapshots at interval={}",
                     uniform_snapshots.size(), sim.transient_output_interval);
        ExportTransientResults(uniform_snapshots);
    }
    else
    {
        // Fixed dt: snapshots are already at output times (old behavior)
        ExportTransientResults(raw_snapshots);
    }
}

// ----------------------------------------------------------------
// Linear interpolation of raw snapshots to a uniform output grid
// ----------------------------------------------------------------
std::vector<TransientSnapshot> TransientSolver::InterpolateSnapshots(
    const std::vector<TransientSnapshot> &raw_snapshots, double t_start, double t_end,
    double output_interval) const
{
    std::vector<TransientSnapshot> result;
    if (raw_snapshots.size() < 2)
    {
        result = raw_snapshots;
        return result;
    }

    // Generate uniform output times
    std::vector<double> output_times;
    for (double t = t_start; t <= t_end + 1e-12 * output_interval; t += output_interval)
    {
        output_times.push_back(t);
    }

    size_t raw_idx = 0;
    for (double t_out : output_times)
    {
        // Find the bracketing interval [raw[j], raw[j+1]] containing t_out
        while (raw_idx + 2 < raw_snapshots.size() &&
               raw_snapshots[raw_idx + 1].time < t_out - 1e-12)
        {
            raw_idx++;
        }

        const auto &s0 = raw_snapshots[raw_idx];
        const auto &s1 = raw_snapshots[raw_idx + 1];
        double dt_seg = s1.time - s0.time;
        double alpha = (dt_seg > 1e-30) ? (t_out - s0.time) / dt_seg : 0.0;
        alpha = std::max(0.0, std::min(1.0, alpha));

        TransientSnapshot snap;
        snap.time = t_out;

        // Helper: lerp two vectors → out = (1-α)*a + α*b
        auto lerp = [&](const mfem::Vector &a, const mfem::Vector &b) -> mfem::Vector
        {
            if (a.Size() == 0)
                return mfem::Vector();
            mfem::Vector out(a.Size());
            for (int i = 0; i < a.Size(); i++)
                out(i) = (1.0 - alpha) * a(i) + alpha * b(i);
            return out;
        };

        snap.voltage = lerp(s0.voltage, s1.voltage);
        snap.temperature = lerp(s0.temperature, s1.temperature);
        snap.displacement = lerp(s0.displacement, s1.displacement);
        result.push_back(std::move(snap));
    }
    return result;
}

// ----------------------------------------------------------------
void TransientSolver::ExportTransientResults(const std::vector<TransientSnapshot> &snapshots)
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

    for (const auto &snap : snapshots)
    {
        times.push_back(snap.time);
        voltages.push_back(snap.voltage);
        temperatures.push_back(snap.temperature);
        displacements.push_back(snap.displacement);
    }

    if (!snapshots.empty())
    {
        exporter.ExportTransient(config_.fe.GetScalarFESpace(), config_.fe.GetVectorFESpace(),
                                 times, voltages, temperatures, displacements);
    }

    logger->info("Transient results exported to: {}", out_dir);
}

}  // namespace fem::transient
