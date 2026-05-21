#include "integrator/PoissonIntegrators.hpp"

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

    // 不加 thread_local，这三个矩阵就是普通的栈上局部变量，每次调用都会重新构造、析构。
    // 加了thread_local 后：
    // 这三个矩阵在每个线程中只构造一次，后续调用会重用它们，避免了频繁的构造和析构开销。
    // 类似于线程内部的全局变量，在一个线程中缓存了这三个矩阵
    // 与多线程安全无关。
    thread_local mfem::DenseMatrix grad_phys_buf;
    thread_local mfem::DenseMatrix dshape_buf;
    thread_local mfem::DenseMatrix inv_jac_buf;
    dshape_buf.SetSize(ndof, dim);
    inv_jac_buf.SetSize(dim);
    grad_phys_buf.SetSize(ndof, dim);

    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);
    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcDShape(ip, dshape_buf);
        mfem::CalcInverse(trans.Jacobian(), inv_jac_buf);
        mfem::Mult(dshape_buf, inv_jac_buf, grad_phys_buf);

        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        // Exploit symmetry: compute upper triangle only
        for (int i = 0; i < ndof; ++i)
        {
            for (int j = i; j < ndof; ++j)
            {
                double grad_dot = 0.0;
                for (int d = 0; d < dim; ++d)
                {
                    grad_dot += grad_phys_buf(i, d) * grad_phys_buf(j, d);
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

    thread_local mfem::Vector shape_buf;
    shape_buf.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape_buf(i);
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

    thread_local mfem::Vector shape_buf;
    shape_buf.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf);
        const double w = ip.weight * trans.Weight() * EvaluateCoefficient(coeff_, trans, ip);

        // Exploit symmetry: compute upper triangle only
        for (int i = 0; i < ndof; ++i)
        {
            for (int j = i; j < ndof; ++j)
            {
                elmat(i, j) += w * shape_buf(i) * shape_buf(j);
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

    thread_local mfem::Vector shape_buf;
    shape_buf.SetSize(ndof);
    const mfem::IntegrationRule &ir = ResolveIntegrationRule(el, trans, IntRule);

    for (int iq = 0; iq < ir.GetNPoints(); ++iq)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf);
        const double rhs_val = EvaluateCoefficient(coeff_, trans, ip);
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            elvect(i) += w * rhs_val * shape_buf(i);
        }
    }
}
}  // namespace fem::integrator
