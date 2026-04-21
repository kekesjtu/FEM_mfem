#include "fem/assembly/ThermalAssembler.hpp"

#include "fem/integrator/PoissonIntegrators.hpp"

namespace fem::assembly
{

ThermalMatrices ThermalAssembler::AssembleMatrices(
    mfem::ParFiniteElementSpace &fespace, mfem::Coefficient &k, mfem::Coefficient &rho_cp,
    std::vector<mfem::ConstantCoefficient> &robin_l_coeffs,
    std::vector<mfem::Array<int>> &robin_markers)
{
    ThermalMatrices m;

    m.K_bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
    m.K_bf->AddDomainIntegrator(new integrator::CustomDiffusionIntegrator(k));
    for (size_t i = 0; i < robin_markers.size(); ++i)
    {
        m.K_bf->AddBoundaryIntegrator(
            new integrator::CustomBoundaryMassIntegrator(robin_l_coeffs[i]), robin_markers[i]);
    }
    m.K_bf->Assemble();
    m.K_bf->Finalize();
    m.K.reset(m.K_bf->ParallelAssemble());

    m.C_bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
    m.C_bf->AddDomainIntegrator(new mfem::MassIntegrator(rho_cp));
    m.C_bf->Assemble();
    m.C_bf->Finalize();
    m.C.reset(m.C_bf->ParallelAssemble());

    return m;
}

mfem::Vector ThermalAssembler::AssembleSourceRHS(
    mfem::ParFiniteElementSpace &fespace, mfem::Coefficient &source,
    std::vector<mfem::ConstantCoefficient> &robin_q_coeffs,
    std::vector<mfem::Array<int>> &robin_markers)
{
    auto lf = std::make_unique<mfem::ParLinearForm>(&fespace);
    lf->AddDomainIntegrator(new integrator::CustomDomainLFIntegrator(source));
    for (size_t i = 0; i < robin_markers.size(); ++i)
    {
        lf->AddBoundaryIntegrator(new integrator::CustomBoundaryLFIntegrator(robin_q_coeffs[i]),
                                  robin_markers[i]);
    }
    lf->Assemble();

    mfem::Vector F(fespace.GetTrueVSize());
    lf->ParallelAssemble(F);
    return F;
}

}  // namespace fem::assembly
