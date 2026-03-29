#pragma once

#include "fem/frontend/Config.hpp"
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

    /// Solve the electrostatic field and return the voltage GridFunction
    void Solve();

    mfem::GridFunction &GetVoltage()
    {
        return voltage_;
    }
    const mfem::GridFunction &GetVoltage() const
    {
        return voltage_;
    }

  private:
    frontend::ProjectConfig &config_;
    mfem::GridFunction voltage_;
    mfem::GridFunction *temperature_gf_ = nullptr;
};

}  // namespace fem::physics
