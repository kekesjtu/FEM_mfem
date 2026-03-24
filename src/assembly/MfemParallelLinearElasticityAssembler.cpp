#include "fem/assembly/MfemParallelLinearElasticityAssembler.hpp"

#include "fem/assembly/LinearElasticityIntegrators.hpp"

#include <stdexcept>

namespace fem::assembly
{
MechanicalSystem MfemParallelLinearElasticityAssembler::Assemble(LinearElasticityInput &input,
                                                                 mfem::GridFunction &initial_guess)
{
#if defined(MFEM_USE_MPI)
    auto *par_space = dynamic_cast<mfem::ParFiniteElementSpace *>(&input.space);
    auto *par_guess = dynamic_cast<mfem::ParGridFunction *>(&initial_guess);
    if (!par_space || !par_guess)
    {
        throw std::runtime_error(
            "MfemParallelLinearElasticityAssembler expects ParFiniteElementSpace/ParGridFunction.");
    }

    MechanicalSystem system;

    system.bilinear = std::make_unique<mfem::ParBilinearForm>(par_space);
    system.bilinear->AddDomainIntegrator(new mfem::ElasticityIntegrator(input.lambda, input.mu));
    for (auto &normal_bc : input.normal_displacement_boundaries)
    {
        if (normal_bc.marker.Size() > 0)
        {
            system.bilinear->AddBoundaryIntegrator(
                new NormalDisplacementPenaltyIntegrator(normal_bc.penalty), normal_bc.marker);
        }
    }
    system.bilinear->Assemble();

    system.linear = std::make_unique<mfem::ParLinearForm>(par_space);
    system.linear->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(input.body_force));
    for (auto &traction_bc : input.traction_boundaries)
    {
        if (traction_bc.marker.Size() > 0)
        {
            system.linear->AddBoundaryIntegrator(
                new mfem::VectorBoundaryLFIntegrator(traction_bc.traction), traction_bc.marker);
        }
    }
    for (auto &pressure_bc : input.pressure_boundaries)
    {
        if (pressure_bc.marker.Size() > 0)
        {
            system.linear->AddBoundaryIntegrator(
                new mfem::VectorBoundaryFluxLFIntegrator(pressure_bc.pressure), pressure_bc.marker);
        }
    }
    if (input.thermal_expansion_secant && input.temperature)
    {
        system.linear->AddDomainIntegrator(
            new ThermalStrainLFIntegrator(input.lambda, input.mu, *input.thermal_expansion_secant,
                                          *input.temperature, input.reference_temperature));
    }
    system.linear->Assemble();

    initial_guess = 0.0;
    initial_guess.ProjectBdrCoefficient(input.dirichlet_displacement, input.essential_bdr_marker);

    par_space->GetEssentialTrueDofs(input.essential_bdr_marker, system.essential_tdofs);
    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
#else
    (void)input;
    (void)initial_guess;
    throw std::runtime_error("MfemParallelLinearElasticityAssembler requires MFEM_USE_MPI.");
#endif
}
}  // namespace fem::assembly
