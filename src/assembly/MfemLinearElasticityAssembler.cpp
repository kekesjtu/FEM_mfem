#include "fem/assembly/MfemLinearElasticityAssembler.hpp"
#include "fem/assembly/LinearElasticityIntegrators.hpp"

namespace fem::assembly
{
MechanicalSystem MfemLinearElasticityAssembler::Assemble(LinearElasticityInput &input,
                                                         mfem::GridFunction &initial_guess)
{
    MechanicalSystem system;

    system.bilinear = std::make_unique<mfem::BilinearForm>(&input.space);
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

    system.linear = std::make_unique<mfem::LinearForm>(&input.space);
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

    input.space.GetEssentialTrueDofs(input.essential_bdr_marker, system.essential_tdofs);
    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
}
}  // namespace fem::assembly
