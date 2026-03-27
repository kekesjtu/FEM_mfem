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

CustomDiffusionIntegrator::CustomDiffusionIntegrator()
    : mfem::DiffusionIntegrator(), coeff_(nullptr)
{
}

CustomDiffusionIntegrator::CustomDiffusionIntegrator(mfem::Coefficient &coeff)
    : mfem::DiffusionIntegrator(coeff), coeff_(&coeff)
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

    mfem::DenseMatrix dshape(ndof, dim);
    mfem::DenseMatrix inv_jac(dim);
    mfem::DenseMatrix grad_phys(ndof, dim);

    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);
    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcDShape(ip, dshape);
        mfem::CalcInverse(trans.Jacobian(), inv_jac);
        mfem::Mult(dshape, inv_jac, grad_phys);

        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        for (int i = 0; i < ndof; ++i)
        {
            for (int j = 0; j < ndof; ++j)
            {
                double grad_dot = 0.0;
                for (int d = 0; d < dim; ++d)
                {
                    grad_dot += grad_phys(i, d) * grad_phys(j, d);
                }
                elmat(i, j) += w * grad_dot;
            }
        }
    }
}

CustomDomainLFIntegrator::CustomDomainLFIntegrator(mfem::Coefficient &coeff)
    : mfem::DomainLFIntegrator(coeff), coeff_(&coeff)
{
}

void CustomDomainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                      mfem::ElementTransformation &trans,
                                                      mfem::Vector &elvect)
{
    const int ndof = el.GetDof();
    elvect.SetSize(ndof);
    elvect = 0.0;

    mfem::Vector shape(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape(i);
        }
    }
}

CustomBoundaryMassIntegrator::CustomBoundaryMassIntegrator(mfem::Coefficient &coeff)
    : mfem::MassIntegrator(coeff), coeff_(&coeff)
{
}

void CustomBoundaryMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                                         mfem::ElementTransformation &trans,
                                                         mfem::DenseMatrix &elmat)
{
    const int ndof = el.GetDof();
    elmat.SetSize(ndof, ndof);
    elmat = 0.0;

    mfem::Vector shape(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape);
        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        for (int i = 0; i < ndof; ++i)
        {
            for (int j = 0; j < ndof; ++j)
            {
                elmat(i, j) += w * shape(i) * shape(j);
            }
        }
    }
}

CustomBoundaryLFIntegrator::CustomBoundaryLFIntegrator(mfem::Coefficient &coeff)
    : mfem::BoundaryLFIntegrator(coeff), coeff_(&coeff)
{
}

void CustomBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                        mfem::ElementTransformation &trans,
                                                        mfem::Vector &elvect)
{
    const int ndof = el.GetDof();
    elvect.SetSize(ndof);
    elvect = 0.0;

    mfem::Vector shape(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape(i);
        }
    }
}
}  // namespace fem::assembly
