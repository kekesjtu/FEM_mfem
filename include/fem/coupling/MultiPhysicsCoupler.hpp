#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/physics/ElectrostaticFieldSolver.hpp"
#include "fem/physics/MechanicalFieldSolver.hpp"
#include "fem/physics/ThermalFieldSolver.hpp"

#include <memory>

namespace fem::coupling
{

/// Multi-physics coupling driver (steady-state coupling strategies).
class MultiPhysicsCoupler
{
  public:
    explicit MultiPhysicsCoupler(frontend::ProjectConfig &config);

    void Solve();

  private:
    void SolveSteadyState();

    void SolveElectroThermalOneWay(std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
                                   std::unique_ptr<physics::ThermalFieldSolver> &t_solver);

    void SolveElectroThermalPicard(std::unique_ptr<physics::ElectrostaticFieldSolver> &e_solver,
                                   std::unique_ptr<physics::ThermalFieldSolver> &t_solver);

    frontend::ProjectConfig &config_;
};

}  // namespace fem::coupling
