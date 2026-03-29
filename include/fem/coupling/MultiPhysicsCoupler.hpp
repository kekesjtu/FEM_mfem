#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"

#include <memory>
#include <vector>

namespace fem::coupling
{

/// Snapshot of all field solutions at one time instant.
struct TransientSnapshot
{
    double time;
    mfem::Vector voltage;
    mfem::Vector temperature;
    mfem::Vector displacement;
};

/// Multi-physics coupling driver.
class MultiPhysicsCoupler
{
  public:
    explicit MultiPhysicsCoupler(frontend::ProjectConfig &config);

    void Solve();

  private:
    void SolveSteadyState();
    void SolveTransient();

    void SolveElectroThermalOneWay(std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
                                   std::unique_ptr<physics::ThermalFieldSolver> &t_solver);

    void SolveElectroThermalPicard(std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
                                   std::unique_ptr<physics::ThermalFieldSolver> &t_solver);

    void ExportTransientResults(const std::vector<TransientSnapshot> &snapshots);

    frontend::ProjectConfig &config_;
};

}  // namespace fem::coupling
