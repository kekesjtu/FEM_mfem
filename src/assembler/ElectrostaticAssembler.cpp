#include "assembler/ElectrostaticAssembler.hpp"

#include "coeff/ElectrostaticCoeffs.hpp"
#include "assembler/LocalAssemblyExecutor.hpp"
#include "integrator/PoissonIntegrators.hpp"

#include <stdexcept>

namespace fem::assembler
{

AssembledSystem ElectrostaticAssembler::Assemble(mfem::ParFiniteElementSpace &space,
                                                 coeff::ElectrostaticCoeffs &coeffs,
                                                 BCbuilder::ScalarBCData &bc,
                                                 mfem::GridFunction &voltage)
{
    auto *sigma = coeffs.GetSigma();
    if (!sigma)
    {
        throw std::runtime_error(
            "ElectrostaticAssembler::Assemble requires an electrical conductivity coefficient");
    }

    AssembledSystem system;

    std::vector<LinearIntegratorFactory> domain_lf;
    domain_lf.push_back([&coeffs]()
                        {
                            return std::make_unique<integrator::CustomDomainLFIntegrator>(
                                coeffs.GetSource());
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
    system.linear = AssembleLinearFormOpenMP(space, domain_lf, boundary_lf);

    std::vector<BilinearIntegratorFactory> domain_bf;
    domain_bf.push_back([sigma]()
                        {
                            return std::make_unique<integrator::CustomDiffusionIntegrator>(
                                *sigma);
                        });

    std::vector<BoundaryBilinearIntegratorSpec> boundary_bf;
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        auto *coeff = &bc.robin_l_coeffs[i];
        auto *marker = &bc.robin_markers[i];
        boundary_bf.push_back(
            {marker,
             [coeff]()
             {
                 return std::make_unique<integrator::CustomBoundaryMassIntegrator>(*coeff);
             }});
    }
    system.bilinear = AssembleBilinearFormOpenMP(space, domain_bf, boundary_bf);
    system.bilinear->Finalize();

    space.GetEssentialTrueDofs(bc.essential_bdr, system.essential_tdofs);
    voltage = 0.0;
    for (size_t i = 0; i < bc.dirichlet_markers.size(); ++i)
    {
        voltage.ProjectBdrCoefficient(bc.dirichlet_coeffs[i], bc.dirichlet_markers[i]);
    }

    system.bilinear->FormLinearSystem(system.essential_tdofs, voltage, *system.linear, system.A,
                                      system.X, system.B);

    return system;
}

}  // namespace fem::assembler
