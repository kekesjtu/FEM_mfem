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

/// Predictor-corrector adaptive time stepper with WRMS error control.
/// Handles backward-Euler time integration, Picard coupling per step,
/// snapshot collection, interpolation, and export.
class TransientSolver
{
  public:
    explicit TransientSolver(frontend::ProjectConfig &config);

    /// Run the full transient simulation.
    void Run();

  private:
    /// Interpolate raw snapshots to a uniform output_interval grid.
    std::vector<TransientSnapshot> InterpolateSnapshots(
        const std::vector<TransientSnapshot> &raw_snapshots, double t_start, double t_end,
        double output_interval) const;

    /// Export transient snapshots (ParaView + txt).
    void ExportTransientResults(const std::vector<TransientSnapshot> &snapshots);

    frontend::ProjectConfig &config_;
};

}  // namespace fem::transient
