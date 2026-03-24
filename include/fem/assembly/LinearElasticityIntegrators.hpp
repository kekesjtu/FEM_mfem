#pragma once

#include "mfem.hpp"

namespace fem::assembly
{
class ThermalStrainLFIntegrator : public mfem::LinearFormIntegrator
{
  public:
    ThermalStrainLFIntegrator(mfem::Coefficient &lambda, mfem::Coefficient &mu,
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

class NormalDisplacementPenaltyIntegrator : public mfem::BilinearFormIntegrator
{
  public:
    explicit NormalDisplacementPenaltyIntegrator(mfem::Coefficient &penalty);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient &penalty_;
};
}  // namespace fem::assembly
