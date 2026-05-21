#pragma once

#include "mfem.hpp"

#include <memory>

namespace fem::frontend
{
struct ProjectConfig;
}

namespace fem::coeff
{

class MechanicalCoeffs
{
  public:
    void BuildCoeffs(frontend::ProjectConfig &config);
    void SetTemperatureField(mfem::GridFunction *temperature_gf);

    mfem::Coefficient *GetLambda() const
    {
        return lambda_.get();
    }

    mfem::Coefficient *GetMu() const
    {
        return mu_.get();
    }

    mfem::VectorCoefficient *GetBodyForce() const
    {
        return body_force_.get();
    }

    mfem::Coefficient *GetAlpha() const
    {
        return alpha_.get();
    }

    mfem::Coefficient *GetTemperature() const
    {
        return temperature_;
    }

    double GetReferenceTemperature() const
    {
        return reference_temperature_value_;
    }

  private:
    std::unique_ptr<mfem::Coefficient> lambda_;
    std::unique_ptr<mfem::Coefficient> mu_;
    std::unique_ptr<mfem::VectorCoefficient> body_force_;
    std::unique_ptr<mfem::Coefficient> alpha_;

    std::unique_ptr<mfem::Coefficient> reference_temperature_coeff_;
    std::unique_ptr<mfem::GridFunctionCoefficient> field_temperature_coeff_;
    
    mfem::Coefficient *temperature_ = nullptr;
    double reference_temperature_value_ = 293.15;
};

}  // namespace fem::coeff
