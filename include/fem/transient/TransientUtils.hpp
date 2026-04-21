#pragma once

#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"
#include "fem/transient/TransientSolver.hpp"
#include "mfem.hpp"

#include <vector>

namespace fem::transient
{

struct FieldSnapshot
{
    double time = 0.0;
    mfem::Vector voltage;
    mfem::Vector temperature;
};

FieldSnapshot SaveFieldSnapshot(double t, physics::ElectrostaticFieldSolver *e_solver,
                                physics::ThermalFieldSolver *t_solver);

void EmitOutputSnapshot(double t_out, const FieldSnapshot &s0, const FieldSnapshot &s1,
                        std::vector<TransientSnapshot> &output_snapshots);

void EmitOutputSnapshot(double t_out, const FieldSnapshot &s0, const FieldSnapshot &s1,
                        const FieldSnapshot &s2, std::vector<TransientSnapshot> &output_snapshots);

struct WRMSResult
{
    double wrms = 0.0;
    double dt_new = 0.0;
    bool accept = true;
};

WRMSResult ComputeAdaptiveError(const mfem::Vector &T_pred, const mfem::Vector &T_corr, double dt,
                                double reltol, double abstol, double eta_safety, double eta_max,
                                double eta_min, double dt_min, int order = 1);

}  // namespace fem::transient
