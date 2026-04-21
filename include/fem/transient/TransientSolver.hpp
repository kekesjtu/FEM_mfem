#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

#include <vector>

namespace fem::transient
{

/// Snapshot of all field solutions at one time instant.
struct TransientSnapshot
{
    double time;
    mfem::Vector voltage;
    mfem::Vector temperature;
    mfem::Vector displacement;
};

/// Predictor-corrector adaptive time stepper with order selection and WRMS error control.
/// Handles backward-Euler / BDF2 time integration, Picard coupling
/// per step, streaming snapshot output, and export.
class TransientSolver
{
  public:
    explicit TransientSolver(frontend::ProjectConfig &config);

    /// Run the full transient simulation.
    void SolveTransient();

  private:
    /// Export transient snapshots (ParaView + txt).
    void ExportTransientResults(std::vector<TransientSnapshot> snapshots);

    frontend::ProjectConfig &config_;
};

}  // namespace fem::transient
