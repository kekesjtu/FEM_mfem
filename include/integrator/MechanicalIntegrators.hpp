#pragma once

#include "mfem.hpp"

namespace fem::integrator
{
class CustomThermalStrainLFIntegrator : public mfem::LinearFormIntegrator
{
  public:
    CustomThermalStrainLFIntegrator(mfem::Coefficient &lambda, mfem::Coefficient &mu,
                                    mfem::Coefficient &alpha_sec, mfem::Coefficient &temperature,
                                    double reference_temperature);

    void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                                mfem::Vector &elvect) override;

  private:
    mfem::Coefficient &lambda_;
    mfem::Coefficient &mu_;
    mfem::Coefficient &alpha_sec_;
    mfem::Coefficient &temperature_;
    double reference_temperature_;
};

class CustomNormalDisplacementPenaltyIntegrator : public mfem::BilinearFormIntegrator
{
  public:
    explicit CustomNormalDisplacementPenaltyIntegrator(mfem::Coefficient &penalty);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient &penalty_;
};
}  // namespace fem::integrator
