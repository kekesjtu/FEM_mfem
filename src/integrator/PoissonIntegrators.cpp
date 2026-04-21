#include "fem/integrator/PoissonIntegrators.hpp"

namespace fem::integrator
{
namespace
{
const mfem::IntegrationRule &ResolveIntegrationRule(const mfem::FiniteElement &el,
                                                    mfem::ElementTransformation &trans,
                                                    const mfem::IntegrationRule *int_rule)
{
    if (int_rule)
    {
        return *int_rule;
    }

    const int order = 2 * el.GetOrder() + trans.OrderW();
    return mfem::IntRules.Get(el.GetGeomType(), order);
}

double EvaluateCoefficient(mfem::Coefficient *coeff, mfem::ElementTransformation &trans,
                           const mfem::IntegrationPoint &ip)
{
    if (!coeff)
    {
        return 1.0;
    }
    return coeff->Eval(trans, ip);
}
}  // namespace

CustomDiffusionIntegrator::CustomDiffusionIntegrator() : coeff_(nullptr)
{
}

CustomDiffusionIntegrator::CustomDiffusionIntegrator(mfem::Coefficient &coeff) : coeff_(&coeff)
{
}

void CustomDiffusionIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                                      mfem::ElementTransformation &trans,
                                                      mfem::DenseMatrix &elmat)
{
    const int ndof = el.GetDof();
    const int dim = el.GetDim();

    elmat.SetSize(ndof, ndof);
    elmat = 0.0;

    dshape_buf_.SetSize(ndof, dim);
    inv_jac_buf_.SetSize(dim);
    grad_phys_buf_.SetSize(ndof, dim);

    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);
    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcDShape(ip, dshape_buf_);
        mfem::CalcInverse(trans.Jacobian(), inv_jac_buf_);
        mfem::Mult(dshape_buf_, inv_jac_buf_, grad_phys_buf_);

        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        // Exploit symmetry: compute upper triangle only
        for (int i = 0; i < ndof; ++i)
        {
            for (int j = i; j < ndof; ++j)
            {
                double grad_dot = 0.0;
                for (int d = 0; d < dim; ++d)
                {
                    grad_dot += grad_phys_buf_(i, d) * grad_phys_buf_(j, d);
                }
                elmat(i, j) += w * grad_dot;
            }
        }
    }

    // Mirror upper triangle to lower
    for (int i = 0; i < ndof; ++i)
        for (int j = 0; j < i; ++j)
            elmat(i, j) = elmat(j, i);
}

CustomDomainLFIntegrator::CustomDomainLFIntegrator(mfem::Coefficient &coeff) : coeff_(&coeff)
{
}

void CustomDomainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                      mfem::ElementTransformation &trans,
                                                      mfem::Vector &elvect)
{
    const int ndof = el.GetDof();
    elvect.SetSize(ndof);
    elvect = 0.0;

    shape_buf_.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf_);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape_buf_(i);
        }
    }
}

CustomBoundaryMassIntegrator::CustomBoundaryMassIntegrator(mfem::Coefficient &coeff)
    : coeff_(&coeff)
{
}

void CustomBoundaryMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                                         mfem::ElementTransformation &trans,
                                                         mfem::DenseMatrix &elmat)
{
    const int ndof = el.GetDof();
    elmat.SetSize(ndof, ndof);
    elmat = 0.0;

    shape_buf_.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf_);
        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        // Exploit symmetry: compute upper triangle only
        for (int i = 0; i < ndof; ++i)
        {
            for (int j = i; j < ndof; ++j)
            {
                elmat(i, j) += w * shape_buf_(i) * shape_buf_(j);
            }
        }
    }

    // Mirror upper triangle to lower
    for (int i = 0; i < ndof; ++i)
        for (int j = 0; j < i; ++j)
            elmat(i, j) = elmat(j, i);
}

CustomBoundaryLFIntegrator::CustomBoundaryLFIntegrator(mfem::Coefficient &coeff) : coeff_(&coeff)
{
}

void CustomBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                        mfem::ElementTransformation &trans,
                                                        mfem::Vector &elvect)
{
    const int ndof = el.GetDof();
    elvect.SetSize(ndof);
    elvect = 0.0;

    shape_buf_.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf_);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape_buf_(i);
        }
    }
}
}  // namespace fem::integrator
