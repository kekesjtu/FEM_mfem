#include "fem/physics/MechanicalFieldSolver.hpp"

#include "fem/assembler/MechanicalAssembler.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

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
        logger->warn("Stiffness matrix not cached, falling back to full Solve()");
        Solve();
        return;
    }

    // Assemble RHS via assembler
    auto lf = assembler::MechanicalAssembler::AssembleRHS(config_.fe.GetVectorFESpace(), coeffs_,
                                                          cached_bc_);

    // Apply essential BCs to initial guess
    displacement_ = 0.0;
    for (size_t i = 0; i < cached_bc_.disp_markers.size(); ++i)
        displacement_.ProjectBdrCoefficient(cached_bc_.disp_coeffs[i], cached_bc_.disp_markers[i]);

    // Form linear system using cached bilinear form
    mfem::OperatorPtr A;
    mfem::Vector X, B;
    cached_bilinear_->FormLinearSystem(cached_bc_.essential_tdofs, displacement_, *lf, A, X, B);

    if (!solver_factorized_)
    {
        cached_solver_ = solver::CreateLinearSolver(config_.simulation.GetSolver("mechanical"),
                                                    1e-12, 1e-12, 5000, 0);
        cached_solver_->Solve(*A, B, X);
        solver_factorized_ = true;
    }
    else
    {
        cached_solver_->Mult(B, X);
    }
    cached_bilinear_->RecoverFEMSolution(X, *lf, displacement_);

    logger->debug("Mechanical: U=[{:.6e}, {:.6e}]", displacement_.Min(), displacement_.Max());
}

}  // namespace fem::physics
