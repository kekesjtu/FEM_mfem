#include "fem/assembler/ElectrostaticAssembler.hpp"

#include "fem/coeff/ElectrostaticCoeffs.hpp"
#include "fem/integrator/PoissonIntegrators.hpp"

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

    system.linear = std::make_unique<mfem::ParLinearForm>(&space);
    system.linear->AddDomainIntegrator(
        new integrator::CustomDomainLFIntegrator(coeffs.GetSource()));
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        system.linear->AddBoundaryIntegrator(
            new integrator::CustomBoundaryLFIntegrator(bc.robin_q_coeffs[i]), bc.robin_markers[i]);
    }
    system.linear->Assemble();

    system.bilinear = std::make_unique<mfem::ParBilinearForm>(&space);
    system.bilinear->AddDomainIntegrator(new integrator::CustomDiffusionIntegrator(*sigma));
    for (size_t i = 0; i < bc.robin_markers.size(); ++i)
    {
        system.bilinear->AddBoundaryIntegrator(
            new integrator::CustomBoundaryMassIntegrator(bc.robin_l_coeffs[i]),
            bc.robin_markers[i]);
    }
    system.bilinear->Assemble();

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