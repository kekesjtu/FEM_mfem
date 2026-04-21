#include "fem/coupling/PicardCoupler.hpp"

#include "fem/log/Logger.hpp"

#include <mpi.h>
#include <cmath>

namespace fem::coupling
{

PicardCoupler::PicardCoupler(const frontend::SimulationConfig &sim_config)
    : max_iterations_(sim_config.picard_max_iterations),
      tolerance_(sim_config.picard_tolerance),
      relaxation_(sim_config.picard_relaxation)
{
}

PicardResult PicardCoupler::RunPicard(physics::ElectrostaticFieldSolver *e_solver,
                                      physics::ThermalFieldSolver *t_solver, bool has_e, bool has_t,
                                      PicardSolveMode mode)
{
    if (has_e && has_t)
    {
        auto logger = fem::log::Get();

        // Lazy-init working vectors on first call (or if DOF count changes)
        const int sz = t_solver->GetTemperature().ParFESpace()->GetTrueVSize();
        if (T_prev_.Size() != sz)
        {
            T_prev_.SetSize(sz);
            delta_.SetSize(sz);
        }

        t_solver->GetTemperature().GetTrueDofs(T_prev_);

        for (int iter = 0; iter < max_iterations_; ++iter)
        {
            e_solver->Solve();
            if (mode == PicardSolveMode::Steady)
                t_solver->SolveSteady();
            else
                t_solver->SolveBDF1();

            // Convergence check (true-dof space)
            t_solver->GetTemperature().GetTrueDofs(delta_);
            delta_ -= T_prev_;
            MPI_Comm comm = t_solver->GetTemperature().ParFESpace()->GetComm();
            double diff_sq = mfem::InnerProduct(comm, delta_, delta_);
            double ref_sq = mfem::InnerProduct(comm, T_prev_, T_prev_);
            double rel_change =
                (ref_sq > 0.0) ? std::sqrt(diff_sq) / std::sqrt(ref_sq) : std::sqrt(diff_sq);

            logger->debug("  Picard iter {}/{}: rel_change={:.6e}", iter + 1, max_iterations_,
                          rel_change);

            if (rel_change < tolerance_)
                return {iter + 1, true};

            // Relaxation
            if (std::abs(relaxation_ - 1.0) > 1e-14)
            {
                mfem::Vector T_cur_true(T_prev_.Size());
                t_solver->GetTemperature().GetTrueDofs(T_cur_true);
                T_cur_true *= relaxation_;
                T_cur_true.Add(1.0 - relaxation_, T_prev_);
                t_solver->GetTemperature().Distribute(T_cur_true);
            }

            t_solver->GetTemperature().GetTrueDofs(T_prev_);
        }
        return {max_iterations_, false};
    }
    else if (has_t)
    {
        if (mode == PicardSolveMode::Steady)
            t_solver->SolveSteady();
        else
            t_solver->SolveBDF1();
        return {1, true};
    }
    else if (has_e)
    {
        e_solver->Solve();
        return {1, true};
    }
    return {0, true};
}

}  // namespace fem::coupling
