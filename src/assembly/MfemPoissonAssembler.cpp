#include "fem/assembly/MfemPoissonAssembler.hpp"

namespace fem::assembly
{
AssembledSystem MfemPoissonAssembler::Assemble(PoissonAssemblyInput &input,
                                               mfem::GridFunction &initial_guess)
{
    AssembledSystem system;

    system.bilinear = std::make_unique<mfem::BilinearForm>(&input.space);
    system.bilinear->AddDomainIntegrator(new mfem::DiffusionIntegrator(input.diffusion));
    if (input.robin_l && input.robin_bdr_marker.Size() > 0)
    {
        system.bilinear->AddBoundaryIntegrator(new mfem::MassIntegrator(*input.robin_l),
                                               input.robin_bdr_marker);
    }
    system.bilinear->Assemble();

    system.linear = std::make_unique<mfem::LinearForm>(&input.space);
    system.linear->AddDomainIntegrator(new mfem::DomainLFIntegrator(input.source));
    if (input.robin_q && input.robin_bdr_marker.Size() > 0)
    {
        system.linear->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*input.robin_q),
                                             input.robin_bdr_marker);
    }
    system.linear->Assemble();

    input.space.GetEssentialTrueDofs(input.essential_bdr_marker, system.essential_tdofs);
    initial_guess = input.dirichlet_value;

    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
}
}  // namespace fem::assembly
