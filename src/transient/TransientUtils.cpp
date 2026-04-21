#include "fem/transient/TransientUtils.hpp"

#include <mpi.h>
#include <algorithm>
#include <cmath>

namespace fem::transient
{

// ----------------------------------------------------------------
// FieldSnapshot helpers
// ----------------------------------------------------------------

FieldSnapshot SaveFieldSnapshot(double t, physics::ElectrostaticFieldSolver *e_solver,
                                physics::ThermalFieldSolver *t_solver)
{
    FieldSnapshot s;
    s.time = t;
    if (e_solver)
        s.voltage = e_solver->GetVoltage();
    if (t_solver)
        s.temperature = t_solver->GetTemperature();
    return s;
}

void EmitOutputSnapshot(double t_out, const FieldSnapshot &s0, const FieldSnapshot &s1,
                        std::vector<TransientSnapshot> &output_snapshots)
{
    double dt_seg = s1.time - s0.time;
    double alpha = (dt_seg > 1e-30) ? (t_out - s0.time) / dt_seg : 0.0;
    alpha = std::max(0.0, std::min(1.0, alpha));

    TransientSnapshot snap;
    snap.time = t_out;
    if (s0.voltage.Size() > 0)
    {
        snap.voltage.SetSize(s0.voltage.Size());
        add(1.0 - alpha, s0.voltage, alpha, s1.voltage, snap.voltage);
    }
    if (s0.temperature.Size() > 0)
    {
        snap.temperature.SetSize(s0.temperature.Size());
        add(1.0 - alpha, s0.temperature, alpha, s1.temperature, snap.temperature);
    }
    output_snapshots.push_back(std::move(snap));
}

// Quadratic (Lagrange) interpolation using 3 snapshots at times t0, t1, t2
void EmitOutputSnapshot(double t_out, const FieldSnapshot &s0, const FieldSnapshot &s1,
                        const FieldSnapshot &s2, std::vector<TransientSnapshot> &output_snapshots)
{
    const double t0 = s0.time, t1 = s1.time, t2 = s2.time;
    // Lagrange basis values at t_out
    const double L0 = ((t_out - t1) * (t_out - t2)) / ((t0 - t1) * (t0 - t2));
    const double L1 = ((t_out - t0) * (t_out - t2)) / ((t1 - t0) * (t1 - t2));
    const double L2 = ((t_out - t0) * (t_out - t1)) / ((t2 - t0) * (t2 - t1));

    TransientSnapshot snap;
    snap.time = t_out;
    if (s0.voltage.Size() > 0)
    {
        snap.voltage.SetSize(s0.voltage.Size());
        // v = L0*v0 + L1*v1 + L2*v2
        add(L0, s0.voltage, L1, s1.voltage, snap.voltage);
        snap.voltage.Add(L2, s2.voltage);
    }
    if (s0.temperature.Size() > 0)
    {
        snap.temperature.SetSize(s0.temperature.Size());
        add(L0, s0.temperature, L1, s1.temperature, snap.temperature);
        snap.temperature.Add(L2, s2.temperature);
    }
    output_snapshots.push_back(std::move(snap));
}

// ----------------------------------------------------------------
// WRMS error from predictor-corrector pair
// ----------------------------------------------------------------

WRMSResult ComputeAdaptiveError(const mfem::Vector &T_pred, const mfem::Vector &T_corr, double dt,
                                double reltol, double abstol, double eta_safety, double eta_max,
                                double eta_min, double dt_min, int order)
{
    int local_N = T_pred.Size();
    double local_sum_sq = 0.0;
    for (int i = 0; i < local_N; i++)
    {
        double e_i = std::abs(T_corr(i) - T_pred(i));
        double weight = abstol + reltol * std::abs(T_corr(i));
        double ratio = e_i / weight;
        local_sum_sq += ratio * ratio;
    }

    double global_sum_sq = local_sum_sq;
    long long global_N = local_N;
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    long long ln = local_N;
    MPI_Allreduce(&ln, &global_N, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    double wrms = std::sqrt(global_sum_sq / global_N);

    double factor = (wrms > 1e-10) ? eta_safety * std::pow(wrms, -1.0 / (order + 1)) : eta_max;
    factor = std::min(factor, eta_max);
    factor = std::max(factor, eta_min);
    double dt_new = dt * factor;

    bool accept = (wrms <= 1.0);
    if (!accept)
        dt_new = std::max(dt_new, dt_min);

    return {wrms, dt_new, accept};
}

}  // namespace fem::transient
