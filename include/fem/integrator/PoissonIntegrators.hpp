#pragma once

#include "mfem.hpp"

namespace fem::integrator
{

//
class CustomDiffusionIntegrator : public mfem::BilinearFormIntegrator
{
  public:
    CustomDiffusionIntegrator();
    explicit CustomDiffusionIntegrator(mfem::Coefficient &coeff);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient *coeff_;
    mutable mfem::DenseMatrix dshape_buf_, inv_jac_buf_, grad_phys_buf_;
};

class CustomDomainLFIntegrator : public mfem::LinearFormIntegrator
{
  public:
    explicit CustomDomainLFIntegrator(mfem::Coefficient &coeff);

    void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                                mfem::Vector &elvect) override;

  private:
    mfem::Coefficient *coeff_;
    mutable mfem::Vector shape_buf_;
};

class CustomBoundaryMassIntegrator : public mfem::BilinearFormIntegrator
{
  public:
    explicit CustomBoundaryMassIntegrator(mfem::Coefficient &coeff);

    void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                               mfem::DenseMatrix &elmat) override;

  private:
    mfem::Coefficient *coeff_;
    mutable mfem::Vector shape_buf_;
};

class CustomBoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
  public:
    explicit CustomBoundaryLFIntegrator(mfem::Coefficient &coeff);

    void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &trans,
                                mfem::Vector &elvect) override;

  private:
    mfem::Coefficient *coeff_;
    mutable mfem::Vector shape_buf_;
};
}  // namespace fem::integrator
