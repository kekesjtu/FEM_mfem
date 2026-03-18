#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/assembly/PoissonIntegrators.hpp"

namespace fem::assembly
{
AssembledSystem MfemPoissonAssembler::Assemble(PoissonAssemblyInput &input,
                                               mfem::GridFunction &initial_guess)
{
    AssembledSystem system;

    system.bilinear = std::make_unique<mfem::BilinearForm>(&input.space);
    system.bilinear->AddDomainIntegrator(new CustomDiffusionIntegrator(input.diffusion));
    for (auto &robin_bc : input.robin_conditions)
    {
        if (robin_bc.marker.Size() > 0)
        {
            system.bilinear->AddBoundaryIntegrator(new CustomBoundaryMassIntegrator(robin_bc.l),
                                                   robin_bc.marker);
        }
    }
    system.bilinear->Assemble();

    system.linear = std::make_unique<mfem::LinearForm>(&input.space);
    system.linear->AddDomainIntegrator(new CustomDomainLFIntegrator(input.source));
    for (auto &robin_bc : input.robin_conditions)
    {
        if (robin_bc.marker.Size() > 0)
        {
            system.linear->AddBoundaryIntegrator(new CustomBoundaryLFIntegrator(robin_bc.q),
                                                 robin_bc.marker);
        }
    }
    system.linear->Assemble();

    input.space.GetEssentialTrueDofs(input.essential_bdr_marker, system.essential_tdofs);
    initial_guess = 0.0;
    if (!input.dirichlet_conditions.empty())
    {
        for (auto &dirichlet_bc : input.dirichlet_conditions)
        {
            if (dirichlet_bc.marker.Size() > 0)
            {
                initial_guess.ProjectBdrCoefficient(dirichlet_bc.value, dirichlet_bc.marker);
            }
        }
    }
    else
    {
        initial_guess = input.dirichlet_value;
    }

    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
}
}  // namespace fem::assembly
