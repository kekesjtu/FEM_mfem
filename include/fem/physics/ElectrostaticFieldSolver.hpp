#pragma once

#include "fem/assembly/ElectrostaticAssembler.hpp"
#include "fem/frontend/Config.hpp"
#include "fem/frontend/Expression.hpp"
#include "fem/solver/ILinearSolver.hpp"
#include "mfem.hpp"

#include <memory>
#include <vector>

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

    struct Coefficients
    {
        std::unique_ptr<mfem::Coefficient> local_sigma;
        mfem::Coefficient *sigma = nullptr;
        mfem::ConstantCoefficient source{0.0};
    };
    Coefficients BuildCoefficients();

    mfem::ParGridFunction &GetVoltage()
    {
        return voltage_;
    }
    const mfem::ParGridFunction &GetVoltage() const
    {
        return voltage_;
    }

  private:
    void BuildBCMarkers();

    frontend::ProjectConfig &config_;
    mfem::ParGridFunction voltage_;
    mfem::GridFunction *temperature_gf_ = nullptr;
    mfem::Coefficient *sigma_ = nullptr;

    /// Parsed expressions for time-varying Dirichlet BCs (index matches dirichlet_bcs).
    std::vector<frontend::Expression> bc_expressions_;

    // --- Cached BC data (markers are immutable; only Dirichlet values update) ---
    assembly::ElectrostaticBCData cached_bc_;
    mfem::Array<int> cached_essential_tdofs_;

    // --- Cached linear solver (reuse across Picard iterations) ---
    std::unique_ptr<solver::ILinearSolver> cached_solver_;
};

}  // namespace fem::physics
