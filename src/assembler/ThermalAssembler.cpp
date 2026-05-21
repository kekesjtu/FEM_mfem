#include "assembler/ThermalAssembler.hpp"

#include "assembler/LocalAssemblyExecutor.hpp"
#include "coeff/ThermalCoeffs.hpp"
#include "integrator/PoissonIntegrators.hpp"

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

    std::vector<BilinearIntegratorFactory> k_domain;
    k_domain.push_back([k]()
                       { return std::make_unique<integrator::CustomDiffusionIntegrator>(*k); });

    std::vector<BoundaryBilinearIntegratorSpec> k_boundary;
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        auto *coeff = &bc.robin_l_coeffs[i];
        auto *marker = &bc.robin_markers[i];
        k_boundary.push_back(
            {marker,
             [coeff]()
             {
                 return std::make_unique<integrator::CustomBoundaryMassIntegrator>(*coeff);
             }});
    }
    m.K_bf = AssembleBilinearFormOpenMP(fespace, k_domain, k_boundary);
    m.K_bf->Finalize();
    m.K.reset(m.K_bf->ParallelAssemble());

    std::vector<BilinearIntegratorFactory> c_domain;
    c_domain.push_back([rho_cp]() { return std::make_unique<mfem::MassIntegrator>(*rho_cp); });
    m.C_bf = AssembleBilinearFormOpenMP(fespace, c_domain, {});
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

    std::vector<LinearIntegratorFactory> domain_lf;
    domain_lf.push_back([source]()
                        {
                            return std::make_unique<integrator::CustomDomainLFIntegrator>(
                                *source);
                        });

    std::vector<BoundaryLinearIntegratorSpec> boundary_lf;
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        auto *coeff = &bc.robin_q_coeffs[i];
        auto *marker = &bc.robin_markers[i];
        boundary_lf.push_back(
            {marker,
             [coeff]()
             {
                 return std::make_unique<integrator::CustomBoundaryLFIntegrator>(*coeff);
             }});
    }
    auto lf = AssembleLinearFormOpenMP(fespace, domain_lf, boundary_lf);

    mfem::Vector F(fespace.GetTrueVSize());
    lf->ParallelAssemble(F);
    return F;
}

}  // namespace fem::assembler
