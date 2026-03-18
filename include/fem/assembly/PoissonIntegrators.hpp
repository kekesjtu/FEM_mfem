#pragma once

#include "mfem.hpp"

namespace fem::assembly
{
class CustomDiffusionIntegrator : public mfem::DiffusionIntegrator
{
  public:
    CustomDiffusionIntegrator();
    explicit CustomDiffusionIntegrator(mfem::Coefficient &coeff);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient *coeff_;
};

class CustomDomainLFIntegrator : public mfem::DomainLFIntegrator
{
  public:
    explicit CustomDomainLFIntegrator(mfem::Coefficient &coeff);

    void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                                mfem::Vector &elvect) override;

  private:
    mfem::Coefficient *coeff_;
};

class CustomBoundaryMassIntegrator : public mfem::MassIntegrator
{
  public:
    explicit CustomBoundaryMassIntegrator(mfem::Coefficient &coeff);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient *coeff_;
};

class CustomBoundaryLFIntegrator : public mfem::BoundaryLFIntegrator
{
  public:
    explicit CustomBoundaryLFIntegrator(mfem::Coefficient &coeff);

    void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                                mfem::Vector &elvect) override;

  private:
    mfem::Coefficient *coeff_;
};
}  // namespace fem::assembly
