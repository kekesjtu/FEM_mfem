#pragma once

#include "fem/BCbuilder/ScalarBCBuilder.hpp"
#include "fem/coeff/ElectrostaticCoeffs.hpp"
#include "fem/frontend/Config.hpp"
#include "fem/solver/ILinearSolver.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::physics
{

class ElectrostaticFieldSolver
{
  public:
    explicit ElectrostaticFieldSolver(frontend::ProjectConfig &config);

    /// Set the temperature GridFunction for coupled sigma(T)
    void SetTemperatureField(mfem::GridFunction *temperature_gf);

    /// Optionally provide an already-built sigma coefficient (e.g. shared with ThermalFieldSolver).
    /// If set, Solve() reuses it instead of constructing a new ExpressionCoefficient internally.
    void SetElectricalConductivity(mfem::Coefficient *sigma);

    /// Update time-varying Dirichlet boundary values for the given time t.
    void UpdateBoundaryConditions(double t);

    void Solve();

    mfem::ParGridFunction &GetVoltage()
    {
        return voltage_;
    }
    const mfem::ParGridFunction &GetVoltage() const
    {
        return voltage_;
    }

    mfem::Coefficient *GetSigma();
    const mfem::Coefficient *GetSigma() const;

  private:
    frontend::ProjectConfig &config_;
    mfem::ParGridFunction voltage_;
    mfem::GridFunction *temperature_gf_ = nullptr;
    mfem::Coefficient *sigma_ = nullptr;

    // --- Field coefficients ---
    coeff::ElectrostaticCoeffs coeffs_;

    // --- Cached BC data (markers are immutable; only Dirichlet values update) ---
    BCbuilder::ScalarBCData cached_bc_;

    // --- Cached linear solver (reuse across Picard iterations) ---
    std::unique_ptr<solver::ILinearSolver> cached_solver_;
};

}  // namespace fem::physics
