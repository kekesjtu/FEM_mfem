#include "fem/post/MechanicalPostProcessor.hpp"

#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

namespace fem::post
{
namespace
{
struct StressState
{
    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    double sxy = 0.0;
    double syz = 0.0;
    double szx = 0.0;
    double von_mises = 0.0;
};

StressState ComputeStressAtPoint(mfem::Mesh &mesh, mfem::GridFunction &displacement,
                                 const LameParameters &lame, int element_id,
                                 const mfem::IntegrationPoint &ip,
                                 const ThermalStrainParameters *thermal)
{
    StressState s;

    auto *tr = mesh.GetElementTransformation(element_id);
    tr->SetIntPoint(&ip);

    mfem::DenseMatrix grad;
    displacement.GetVectorGradient(*tr, grad);

    const double lambda = lame.lambda.Eval(*tr, ip);
    const double mu = lame.mu.Eval(*tr, ip);

    const int dim = mesh.Dimension();
    const double exx = grad(0, 0);
    const double eyy = (dim > 1) ? grad(1, 1) : 0.0;
    const double ezz = (dim > 2) ? grad(2, 2) : 0.0;
    const double exy = (dim > 1) ? 0.5 * (grad(0, 1) + grad(1, 0)) : 0.0;
    const double eyz = (dim > 2) ? 0.5 * (grad(1, 2) + grad(2, 1)) : 0.0;
    const double ezx = (dim > 2) ? 0.5 * (grad(2, 0) + grad(0, 2)) : 0.0;

    const double trace_eps = exx + eyy + ezz;

    s.sxx = lambda * trace_eps + 2.0 * mu * exx;
    s.syy = lambda * trace_eps + 2.0 * mu * eyy;
    s.szz = lambda * trace_eps + 2.0 * mu * ezz;
    s.sxy = 2.0 * mu * exy;
    s.syz = 2.0 * mu * eyz;
    s.szx = 2.0 * mu * ezx;

    if (thermal)
    {
        const double alpha_sec = thermal->alpha_sec.Eval(*tr, ip);
        const double temperature = thermal->temperature.Eval(*tr, ip);
        const double delta_t = temperature - thermal->reference_temperature;
        const double thermal_shift = (3.0 * lambda + 2.0 * mu) * alpha_sec * delta_t;

        s.sxx -= thermal_shift;
        s.syy -= thermal_shift;
        s.szz -= thermal_shift;
    }

    const double vm_sq =
        0.5 * ((s.sxx - s.syy) * (s.sxx - s.syy) + (s.syy - s.szz) * (s.syy - s.szz) +
               (s.szz - s.sxx) * (s.szz - s.sxx)) +
        3.0 * (s.sxy * s.sxy + s.syz * s.syz + s.szx * s.szx);
    s.von_mises = std::sqrt(std::max(0.0, vm_sq));

    return s;
}

StressState ComputeElementStress(mfem::Mesh &mesh, mfem::GridFunction &displacement,
                                 const LameParameters &lame, int element_id,
                                 const ThermalStrainParameters *thermal)
{
    const auto &center = mfem::Geometries.GetCenter(mesh.GetElementBaseGeometry(element_id));
    return ComputeStressAtPoint(mesh, displacement, lame, element_id, center, thermal);
}
}  // namespace

void MechanicalPostProcessor::FillVonMisesNodalField(mfem::Mesh &mesh,
                                                     mfem::GridFunction &displacement,
                                                     const LameParameters &lame,
                                                     mfem::FiniteElementSpace &scalar_space,
                                                     mfem::GridFunction &von_mises,
                                                     const ThermalStrainParameters *thermal)
{
    von_mises = 0.0;

    std::vector<double> weighted_sum(scalar_space.GetNDofs(), 0.0);
    std::vector<double> shape_sum(scalar_space.GetNDofs(), 0.0);

    mfem::Array<int> vdofs;
    for (int e = 0; e < mesh.GetNE(); ++e)
    {
        const auto *fe = scalar_space.GetFE(e);
        auto *tr = mesh.GetElementTransformation(e);
        if (!fe || !tr)
        {
            continue;
        }

        scalar_space.GetElementVDofs(e, vdofs);
        mfem::Vector shape(fe->GetDof());

        const int order = 2 * fe->GetOrder() + tr->OrderW();
        const mfem::IntegrationRule &ir = mfem::IntRules.Get(fe->GetGeomType(), order);

        for (int q = 0; q < ir.GetNPoints(); ++q)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);
            const auto stress = ComputeStressAtPoint(mesh, displacement, lame, e, ip, thermal);

            tr->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);
            const double w = ip.weight * std::abs(tr->Weight());

            for (int i = 0; i < vdofs.Size(); ++i)
            {
                int dof = vdofs[i];
                if (dof < 0)
                {
                    dof = -1 - dof;
                }
                if (dof < 0 || dof >= scalar_space.GetNDofs())
                {
                    continue;
                }

                const double phi = shape(i);
                weighted_sum[dof] += w * phi * stress.von_mises;
                shape_sum[dof] += w * phi;
            }
        }
    }

    const double eps = 1.0e-30;
    for (int d = 0; d < scalar_space.GetNDofs(); ++d)
    {
        if (shape_sum[d] <= eps)
        {
            continue;
        }
        von_mises(d) = weighted_sum[d] / shape_sum[d];
    }
}

void MechanicalPostProcessor::ExportElementStressCsv(const std::string &csv_path, mfem::Mesh &mesh,
                                                     mfem::GridFunction &displacement,
                                                     const LameParameters &lame,
                                                     const ThermalStrainParameters *thermal)
{
    std::filesystem::create_directories(std::filesystem::path(csv_path).parent_path());

    std::ofstream out(csv_path);
    out << "element_id,attribute,cx,cy,cz,sxx,syy,szz,sxy,syz,szx,von_mises\n";

    for (int e = 0; e < mesh.GetNE(); ++e)
    {
        auto *tr = mesh.GetElementTransformation(e);
        const auto &center = mfem::Geometries.GetCenter(mesh.GetElementBaseGeometry(e));
        tr->SetIntPoint(&center);

        mfem::Vector c;
        tr->Transform(center, c);
        const auto stress = ComputeElementStress(mesh, displacement, lame, e, thermal);
        const int attr = mesh.GetAttribute(e);

        const double cx = (c.Size() > 0) ? c(0) : 0.0;
        const double cy = (c.Size() > 1) ? c(1) : 0.0;
        const double cz = (c.Size() > 2) ? c(2) : 0.0;

        out << e << ',' << attr << ',' << cx << ',' << cy << ',' << cz << ',' << stress.sxx << ','
            << stress.syy << ',' << stress.szz << ',' << stress.sxy << ',' << stress.syz << ','
            << stress.szx << ',' << stress.von_mises << '\n';
    }
}
}  // namespace fem::post
