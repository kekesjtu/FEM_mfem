#pragma once

#include "fem/BCbuilder/MechanicalBCBuilder.hpp"
#include "fem/coeff/MechanicalCoeffs.hpp"
#include "fem/frontend/Config.hpp"
#include "fem/solver/ILinearSolver.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::physics
{

class MechanicalFieldSolver
{
  public:
    explicit MechanicalFieldSolver(frontend::ProjectConfig &config);

    /// Optionally set an external temperature GridFunction for thermal expansion coupling.
    /// If not set, the solver uses the reference temperature everywhere, so thermal strain is zero.
    void SetTemperatureField(mfem::GridFunction *temperature_gf);

    /// Full solve (assemble + solve). Used for one-off or steady-state.
    void Solve();

    /// Pre-assemble and cache the stiffness matrix (LHS).
    /// Call once before a batch of SolveRHSOnly() calls.
    void CacheStiffnessMatrix();

    /// Solve with cached stiffness matrix; only reassembles the RHS (thermal strain).
    /// CacheStiffnessMatrix() must have been called first.
    void SolveRHSOnly();

    mfem::ParGridFunction &GetDisplacement()
    {
        return displacement_;
    }
    const mfem::ParGridFunction &GetDisplacement() const
    {
        return displacement_;
    }

  private:
    frontend::ProjectConfig &config_;
    mfem::ParGridFunction displacement_;

    // --- Constant coefficients (built once) ---
    coeff::MechanicalCoeffs coeffs_;

    // --- Cached BC data (built once) ---
    BCbuilder::MechanicalBCData cached_bc_;

    // --- Cached stiffness matrix data ---
    bool stiffness_cached_ = false;
    std::unique_ptr<mfem::ParBilinearForm> cached_bilinear_;

    // Cached linear solver (reuse factorization across snapshots)
    std::unique_ptr<solver::ILinearSolver> cached_solver_;
    bool solver_factorized_ = false;
};

}  // namespace fem::physics
