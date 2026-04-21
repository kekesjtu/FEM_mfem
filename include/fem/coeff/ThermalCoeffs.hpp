#pragma once

#include "fem/coeff/CoefficientManager.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::frontend
{
struct ProjectConfig;
}

namespace fem::coeff
{

class ThermalCoeffs
{
  public:
    void BuildCoeffs(frontend::ProjectConfig &config, mfem::GridFunction * const *voltage_gf,
                     mfem::Coefficient * const *sigma);

    mfem::Coefficient *GetK() const
    {
        return k_.get();
    }

    mfem::Coefficient *GetRhoCp() const
    {
        return rho_cp_.get();
    }

    mfem::Coefficient *GetSource() const
    {
        return total_source_.get();
    }

  private:
    std::unique_ptr<mfem::Coefficient> k_;
    std::unique_ptr<mfem::Coefficient> rho_cp_;
    std::unique_ptr<mfem::Coefficient> base_source_;
    std::unique_ptr<coeff::JouleHeatingCoefficient> joule_source_;
    std::unique_ptr<mfem::Coefficient> total_source_;
    mfem::GridFunction * const *voltage_gf_ = nullptr;
    mfem::Coefficient * const *sigma_ = nullptr;
};

}  // namespace fem::coeff
