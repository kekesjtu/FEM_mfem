#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::physics
{

class MechanicalFieldSolver
{
  public:
    explicit MechanicalFieldSolver(frontend::ProjectConfig &config);

    /// Set temperature GridFunction for thermal expansion coupling
    void SetTemperatureField(mfem::GridFunction *temperature_gf);

    void Solve();

    mfem::GridFunction &GetDisplacement()
    {
        return displacement_;
    }
    const mfem::GridFunction &GetDisplacement() const
    {
        return displacement_;
    }

  private:
    frontend::ProjectConfig &config_;
    mfem::GridFunction displacement_;
    mfem::GridFunction *temperature_gf_ = nullptr;
};

}  // namespace fem::physics
