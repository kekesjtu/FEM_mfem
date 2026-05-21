#include "integrator/MechanicalIntegrators.hpp"

#include <cmath>

namespace fem::integrator
{
CustomThermalStrainLFIntegrator::CustomThermalStrainLFIntegrator(mfem::Coefficient &lambda,
                                                                 mfem::Coefficient &mu,
                                                                 mfem::Coefficient &alpha_sec,
                                                                 mfem::Coefficient &temperature,
                                                                 double reference_temperature)
    : lambda_(lambda),
      mu_(mu),
      alpha_sec_(alpha_sec),
      temperature_(temperature),
      reference_temperature_(reference_temperature)
{
}

void CustomThermalStrainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                             mfem::ElementTransformation &trans,
                                                             mfem::Vector &elvect)
{
    const int dim = el.GetDim();
    const int ndof = el.GetDof();
    elvect.SetSize(dim * ndof);
    elvect = 0.0;

    // 不加 thread_local，这三个矩阵就是普通的栈上局部变量，每次调用都会重新构造、析构。
    // 加了thread_local 后：
    // 这三个矩阵在每个线程中只构造一次，后续调用会重用它们，避免了频繁的构造和析构开销。
    // 相当于在一个线程中缓存了这三个矩阵
    // 与多线程安全无关。
    thread_local mfem::DenseMatrix dshape_buf;
    thread_local mfem::DenseMatrix inv_jac_buf;
    thread_local mfem::DenseMatrix grad_phys_buf;
    dshape_buf.SetSize(ndof, dim);
    inv_jac_buf.SetSize(dim);
    grad_phys_buf.SetSize(ndof, dim);

    const mfem::IntegrationRule *ir = IntRule;
    if (!ir)
    {
        const int order = 2 * el.GetOrder() + trans.OrderW();
        ir = &mfem::IntRules.Get(el.GetGeomType(), order);
    }

    for (int q = 0; q < ir->GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir->IntPoint(q);
        trans.SetIntPoint(&ip);

        el.CalcDShape(ip, dshape_buf);
        mfem::CalcInverse(trans.Jacobian(), inv_jac_buf);
        mfem::Mult(dshape_buf, inv_jac_buf, grad_phys_buf);

        const double lambda = lambda_.Eval(trans, ip);
        const double mu = mu_.Eval(trans, ip);
        const double alpha_sec = alpha_sec_.Eval(trans, ip);
        const double delta_t = temperature_.Eval(trans, ip) - reference_temperature_;
        const double epsilon_th = alpha_sec * delta_t;
        const double beta = (3.0 * lambda + 2.0 * mu) * epsilon_th;
        const double w = ip.weight * trans.Weight();

        for (int i = 0; i < ndof; ++i)
        {
            for (int d = 0; d < dim; ++d)
            {
                const int vdof = i + d * ndof;
                elvect(vdof) += w * beta * grad_phys_buf(i, d);
            }
        }
    }
}

CustomNormalDisplacementPenaltyIntegrator::CustomNormalDisplacementPenaltyIntegrator(
    mfem::Coefficient &penalty)
    : penalty_(penalty)
{
}

void CustomNormalDisplacementPenaltyIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el, mfem::ElementTransformation &trans, mfem::DenseMatrix &elmat)
{
    const int ndof = el.GetDof();
    const int space_dim = trans.GetSpaceDim();
    elmat.SetSize(space_dim * ndof, space_dim * ndof);
    elmat = 0.0;

    thread_local mfem::Vector shape_buf;
    thread_local mfem::Vector normal_buf;
    shape_buf.SetSize(ndof);
    normal_buf.SetSize(space_dim);

    const mfem::IntegrationRule *ir = IntRule;
    if (!ir)
    {
        const int order = 2 * el.GetOrder() + trans.OrderW();
        ir = &mfem::IntRules.Get(el.GetGeomType(), order);
    }

    for (int q = 0; q < ir->GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir->IntPoint(q);
        trans.SetIntPoint(&ip);

        el.CalcShape(ip, shape_buf);

        mfem::CalcOrtho(trans.Jacobian(), normal_buf);
        const double normal_norm = normal_buf.Norml2();
        if (normal_norm <= 0.0)
        {
            continue;
        }
        normal_buf /= normal_norm;

        const double penalty = penalty_.Eval(trans, ip);
        const double w = ip.weight * std::abs(trans.Weight()) * penalty;

        for (int i = 0; i < ndof; ++i)
        {
            for (int j = 0; j < ndof; ++j)
            {
                const double s = w * shape_buf(i) * shape_buf(j);
                for (int d = 0; d < space_dim; ++d)
                {
                    const int row = i + d * ndof;
                    for (int e = 0; e < space_dim; ++e)
                    {
                        const int col = j + e * ndof;
                        elmat(row, col) += s * normal_buf(d) * normal_buf(e);
                    }
                }
            }
        }
    }
}
}  // namespace fem::integrator
