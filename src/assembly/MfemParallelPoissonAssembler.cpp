#include "fem/assembly/MfemParallelPoissonAssembler.hpp"

#include "fem/assembly/PoissonIntegrators.hpp"

#include <stdexcept>

namespace fem::assembly
{
AssembledSystem MfemParallelPoissonAssembler::Assemble(PoissonAssemblyInput &input,
                                                       mfem::GridFunction &initial_guess)
{
#if defined(MFEM_USE_MPI)
    auto *par_space = dynamic_cast<mfem::ParFiniteElementSpace *>(&input.space);
    auto *par_guess = dynamic_cast<mfem::ParGridFunction *>(&initial_guess);
    if (!par_space || !par_guess)
    {
        throw std::runtime_error(
            "MfemParallelPoissonAssembler expects ParFiniteElementSpace/ParGridFunction.");
    }

    AssembledSystem system;

    system.bilinear = std::make_unique<mfem::ParBilinearForm>(par_space);
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

    system.linear = std::make_unique<mfem::ParLinearForm>(par_space);
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

    par_space->GetEssentialTrueDofs(input.essential_bdr_marker, system.essential_tdofs);
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
#else
    (void)input;
    (void)initial_guess;
    throw std::runtime_error("MfemParallelPoissonAssembler requires MFEM_USE_MPI.");
#endif
}
}  // namespace fem::assembly
