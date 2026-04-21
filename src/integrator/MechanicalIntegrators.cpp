#include "fem/integrator/MechanicalIntegrators.hpp"

#include <cmath>

namespace fem::integrator
{
ThermalStrainLFIntegrator::ThermalStrainLFIntegrator(mfem::Coefficient &lambda,
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

void ThermalStrainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                       mfem::ElementTransformation &trans,
                                                       mfem::Vector &elvect)
{
    const int dim = el.GetDim();
    const int ndof = el.GetDof();
    elvect.SetSize(dim * ndof);
    elvect = 0.0;

    dshape_buf_.SetSize(ndof, dim);
    inv_jac_buf_.SetSize(dim);
    grad_phys_buf_.SetSize(ndof, dim);

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

        el.CalcDShape(ip, dshape_buf_);
        mfem::CalcInverse(trans.Jacobian(), inv_jac_buf_);
        mfem::Mult(dshape_buf_, inv_jac_buf_, grad_phys_buf_);

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
                elvect(vdof) += w * beta * grad_phys_buf_(i, d);
            }
        }
    }
}

NormalDisplacementPenaltyIntegrator::NormalDisplacementPenaltyIntegrator(mfem::Coefficient &penalty)
    : penalty_(penalty)
{
}

void NormalDisplacementPenaltyIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                                                mfem::ElementTransformation &trans,
                                                                mfem::DenseMatrix &elmat)
{
    const int ndof = el.GetDof();
    const int space_dim = trans.GetSpaceDim();
    elmat.SetSize(space_dim * ndof, space_dim * ndof);
    elmat = 0.0;

    shape_buf_.SetSize(ndof);
    normal_buf_.SetSize(space_dim);

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

        el.CalcShape(ip, shape_buf_);

        mfem::CalcOrtho(trans.Jacobian(), normal_buf_);
        const double normal_norm = normal_buf_.Norml2();
        if (normal_norm <= 0.0)
        {
            continue;
        }
        normal_buf_ /= normal_norm;

        const double penalty = penalty_.Eval(trans, ip);
        const double w = ip.weight * std::abs(trans.Weight()) * penalty;

        for (int i = 0; i < ndof; ++i)
        {
            for (int j = 0; j < ndof; ++j)
            {
                const double s = w * shape_buf_(i) * shape_buf_(j);
                for (int d = 0; d < space_dim; ++d)
                {
                    const int row = i + d * ndof;
                    for (int e = 0; e < space_dim; ++e)
                    {
                        const int col = j + e * ndof;
                        elmat(row, col) += s * normal_buf_(d) * normal_buf_(e);
                    }
                }
            }
        }
    }
}
}  // namespace fem::integrator
