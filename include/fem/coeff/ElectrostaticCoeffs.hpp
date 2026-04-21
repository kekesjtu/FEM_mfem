#pragma once

#include "fem/BCbuilder/ScalarBCBuilder.hpp"
#include "fem/frontend/Expression.hpp"
#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::frontend
{
struct ProjectConfig;
}

namespace fem::coeff
{

class ElectrostaticCoeffs
{
  public:
    void BuildCoeffs(frontend::ProjectConfig &config, mfem::GridFunction *const *temperature_gf,
                     mfem::Coefficient *const *sigma);
    void UpdateBoundaryConditions(const frontend::ScalarFieldConfig &field_config, double t,
                                  BCbuilder::ScalarBCData &bc) const;

    mfem::Coefficient *GetSigma() const
    {
        if (sigma_ && *sigma_)
            return *sigma_;
        return owned_sigma_.get();
    }

    mfem::ConstantCoefficient &GetSource()
    {
        return source_;
    }

  private:
    std::unique_ptr<mfem::Coefficient> owned_sigma_;
    mfem::Coefficient *const *sigma_ = nullptr;
    mfem::GridFunction *const *temperature_gf_ = nullptr;
    mfem::ConstantCoefficient source_{0.0};
    std::vector<frontend::Expression> bc_expressions_;
};

}  // namespace fem::coeff
