#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::physics
{

class ThermalFieldSolver
{
  public:
    explicit ThermalFieldSolver(frontend::ProjectConfig &config);

    /// Set voltage GridFunction for Joule heating source Q = sigma * |grad V|^2
    void SetVoltageField(mfem::GridFunction *voltage_gf);

    /// Set sigma coefficient for Joule heating (needed for coupled problems)
    void SetElectricalConductivity(mfem::Coefficient *sigma);

    /// Enable backward-Euler transient mode for subsequent Solve() calls.
    void EnableTransient(double dt, mfem::GridFunction *T_old);

    void Solve();

    mfem::GridFunction &GetTemperature()
    {
        return temperature_;
    }
    const mfem::GridFunction &GetTemperature() const
    {
        return temperature_;
    }

  private:
    frontend::ProjectConfig &config_;
    mfem::GridFunction temperature_;
    mfem::GridFunction *voltage_gf_ = nullptr;
    mfem::Coefficient *sigma_ = nullptr;
    double transient_dt_ = 0.0;
    mfem::GridFunction *T_old_ = nullptr;
};

}  // namespace fem::physics
