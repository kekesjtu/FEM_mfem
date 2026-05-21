#pragma once

#include "frontend/Config.hpp"
#include "physics/ElectrostaticFieldSolver.hpp"
#include "physics/ThermalFieldSolver.hpp"
#include "mfem.hpp"

namespace fem::coupling
{

struct PicardResult
{
    int iterations = 0;
    bool converged = false;
};

enum class PicardSolveMode
{
    Steady,
    BackwardEuler
};

class PicardCoupler
{
  public:
    explicit PicardCoupler(const frontend::SimulationConfig &sim_config);

    PicardResult RunPicard(physics::ElectrostaticFieldSolver *e_solver,
                           physics::ThermalFieldSolver *t_solver, bool has_e, bool has_t,
                           PicardSolveMode mode);

  private:
    int max_iterations_;
    double tolerance_;
    double relaxation_;
    mfem::Vector T_prev_true_;
    mfem::Vector delta_;
};

}  // namespace fem::coupling
