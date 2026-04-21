#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/solver/ILinearSolver.hpp"
#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::physics
{

class MechanicalFieldSolver
{
  public:
    explicit MechanicalFieldSolver(frontend::ProjectConfig &config);

    /// Set temperature GridFunction for thermal expansion coupling
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
    void BuildBCMarkers();
    std::unique_ptr<mfem::Coefficient> BuildTempCoefficient();

    frontend::ProjectConfig &config_;
    mfem::ParGridFunction displacement_;
    mfem::GridFunction *temperature_gf_ = nullptr;

    // --- Constant coefficients (built once) ---
    std::unique_ptr<mfem::Coefficient> cached_lambda_;
    std::unique_ptr<mfem::Coefficient> cached_mu_;
    std::unique_ptr<mfem::VectorCoefficient> cached_body_force_;
    std::unique_ptr<mfem::Coefficient> cached_alpha_;

    // --- Cached BC markers (built once) ---
    mfem::Array<int> cached_essential_bdr_;
    mfem::Array<int> cached_essential_tdofs_;
    std::vector<mfem::VectorConstantCoefficient> cached_disp_coeffs_;
    std::vector<mfem::Array<int>> cached_disp_markers_;
    std::vector<mfem::ConstantCoefficient> cached_penalty_coeffs_;
    std::vector<mfem::Array<int>> cached_nd_markers_;
    std::vector<mfem::ConstantCoefficient> cached_pressure_coeffs_;
    std::vector<mfem::Array<int>> cached_pressure_markers_;

    // --- Cached stiffness matrix data ---
    bool stiffness_cached_ = false;
    std::unique_ptr<mfem::ParBilinearForm> cached_bilinear_;

    // Cached linear solver (reuse factorization across snapshots)
    std::unique_ptr<solver::ILinearSolver> cached_solver_;
    bool solver_factorized_ = false;
};

}  // namespace fem::physics
