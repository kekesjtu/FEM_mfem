#include "physics/MechanicalFieldSolver.hpp"

#include "assembler/MechanicalAssembler.hpp"
#include "log/Logger.hpp"
#include "solver/LinearSolverFactory.hpp"

#include <stdexcept>

namespace fem::physics
{

MechanicalFieldSolver::MechanicalFieldSolver(frontend::ProjectConfig &config)
    : config_(config), displacement_(&config.fe.GetVectorFESpace())
{
    displacement_ = 0.0;
    coeffs_.BuildCoeffs(config_);

    // Build and cache BC markers once
    cached_bc_ = BCbuilder::MechanicalBCBuilder::BuildBC(config_);
}

void MechanicalFieldSolver::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    coeffs_.SetTemperatureField(temperature_gf);
}

void MechanicalFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->debug("Mechanical full solve");

    // Use CacheStiffnessMatrix + SolveRHSOnly pattern to avoid redundant assembly/factorization
    CacheStiffnessMatrix();
    SolveRHSOnly();
}

// ------------------------------------------------------------------
// CacheStiffnessMatrix — delegate to MechanicalAssembler
// ------------------------------------------------------------------
void MechanicalFieldSolver::CacheStiffnessMatrix()
{
    if (stiffness_cached_)
        return;

    auto logger = fem::log::Get();
    logger->debug("Caching mechanical stiffness matrix...");

    cached_bilinear_ = assembler::MechanicalAssembler::AssembleStiffness(
        config_.fe.GetVectorFESpace(), coeffs_, cached_bc_);

    stiffness_cached_ = true;
    solver_factorized_ = false;
    cached_solver_.reset();
    logger->debug("Stiffness matrix cached.");
}

// ------------------------------------------------------------------
// SolveRHSOnly — delegate RHS assembly to MechanicalAssembler, solve with cached K
// ------------------------------------------------------------------
void MechanicalFieldSolver::SolveRHSOnly()
{
    auto logger = fem::log::Get();
    logger->debug("Mechanical solve (cached K, RHS-only)");

    if (!stiffness_cached_)
    {
        throw std::runtime_error(
            "MechanicalFieldSolver::SolveRHSOnly requires CacheStiffnessMatrix() first.");
    }

    // Assemble RHS via assembler
    auto lf = assembler::MechanicalAssembler::AssembleRHS(config_.fe.GetVectorFESpace(), coeffs_,
                                                          cached_bc_);

    // The first solve starts from zero; subsequent RHS-only solves reuse the previous
    // displacement as an initial guess for the iterative solver.
    if (!solver_factorized_)
    {
        displacement_ = 0.0;
    }
    for (size_t i = 0; i < cached_bc_.disp_markers.size(); ++i)
    {
        displacement_.ProjectBdrCoefficient(cached_bc_.disp_coeffs[i], cached_bc_.disp_markers[i]);
    }

    // Form linear system using cached bilinear form
    mfem::OperatorPtr A;
    mfem::Vector X, B;
    cached_bilinear_->FormLinearSystem(cached_bc_.essential_tdofs, displacement_, *lf, A, X, B);

    if (!solver_factorized_)
    {
        logger->info("Mechanical linear solve starts: solver='{}', true_ndof={}",
                     config_.simulation.GetSolver("mechanical"),
                     config_.fe.GetVectorFESpace().GlobalTrueVSize());
        cached_solver_ =
            solver::CreateLinearSolver(config_.simulation.GetSolver("mechanical"), 1e-8, 1e-12,
                                       2000, 1, &config_.fe.GetVectorFESpace());
        cached_solver_->Solve(*A, B, X);
        solver_factorized_ = true;
    }
    else
    {
        logger->info("Mechanical linear solve starts: solver='{}', cached operator/RHS",
                     config_.simulation.GetSolver("mechanical"));
        cached_solver_->Mult(B, X);
    }
    cached_bilinear_->RecoverFEMSolution(X, *lf, displacement_);

    logger->info("Mechanical linear solve done.");
    logger->debug("Mechanical: U=[{:.6e}, {:.6e}]", displacement_.Min(), displacement_.Max());
}

}  // namespace fem::physics
