#include "fem/assembler/ThermalAssembler.hpp"

#include "fem/coeff/ThermalCoeffs.hpp"
#include "fem/integrator/PoissonIntegrators.hpp"

#include <stdexcept>

namespace fem::assembler
{

ThermalMatrices ThermalAssembler::AssembleMatrices(mfem::ParFiniteElementSpace &fespace,
                                                   const coeff::ThermalCoeffs &coeffs,
                                                   BCbuilder::ScalarBCData &bc)
{
    auto *k = coeffs.GetK();
    auto *rho_cp = coeffs.GetRhoCp();
    if (!k || !rho_cp)
    {
        throw std::runtime_error(
            "ThermalAssembler::AssembleMatrices requires k and rho_cp coefficients");
    }

    ThermalMatrices m;

    m.K_bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
    m.K_bf->AddDomainIntegrator(new integrator::CustomDiffusionIntegrator(*k));
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        m.K_bf->AddBoundaryIntegrator(
            new integrator::CustomBoundaryMassIntegrator(bc.robin_l_coeffs[i]),
            bc.robin_markers[i]);
    }
    m.K_bf->Assemble();
    m.K_bf->Finalize();
    m.K.reset(m.K_bf->ParallelAssemble());

    m.C_bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
    m.C_bf->AddDomainIntegrator(new mfem::MassIntegrator(*rho_cp));
    m.C_bf->Assemble();
    m.C_bf->Finalize();
    m.C.reset(m.C_bf->ParallelAssemble());

    return m;
}

mfem::Vector ThermalAssembler::AssembleSourceRHS(mfem::ParFiniteElementSpace &fespace,
                                                 const coeff::ThermalCoeffs &coeffs,
                                                 BCbuilder::ScalarBCData &bc)
{
    auto *source = coeffs.GetSource();
    if (!source)
    {
        throw std::runtime_error(
            "ThermalAssembler::AssembleSourceRHS requires a heat source coefficient");
    }

    auto lf = std::make_unique<mfem::ParLinearForm>(&fespace);
    lf->AddDomainIntegrator(new integrator::CustomDomainLFIntegrator(*source));
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        lf->AddBoundaryIntegrator(new integrator::CustomBoundaryLFIntegrator(bc.robin_q_coeffs[i]),
                                  bc.robin_markers[i]);
    }
    lf->Assemble();

    mfem::Vector F(fespace.GetTrueVSize());
    lf->ParallelAssemble(F);
    return F;
}

}  // namespace fem::assembler
