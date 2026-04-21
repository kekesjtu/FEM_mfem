#include "fem/assembler/MechanicalAssembler.hpp"

#include "fem/coeff/MechanicalCoeffs.hpp"
#include "fem/integrator/MechanicalIntegrators.hpp"

#include <stdexcept>

namespace fem::assembler
{

std::unique_ptr<mfem::ParBilinearForm> MechanicalAssembler::AssembleStiffness(
    mfem::ParFiniteElementSpace &space, const coeff::MechanicalCoeffs &coeffs,
    BCbuilder::MechanicalBCData &bc)
{
    auto *lambda = coeffs.GetLambda();
    auto *mu = coeffs.GetMu();
    if (!lambda || !mu)
    {
        throw std::runtime_error(
            "MechanicalAssembler::AssembleStiffness requires lambda and mu coefficients");
    }

    auto bf = std::make_unique<mfem::ParBilinearForm>(&space);
    bf->AddDomainIntegrator(new mfem::ElasticityIntegrator(*lambda, *mu));
    for (size_t i = 0; i < bc.nd_markers.size(); ++i)
    {
        bf->AddBoundaryIntegrator(
            new integrator::NormalDisplacementPenaltyIntegrator(bc.penalty_coeffs[i]),
            bc.nd_markers[i]);
    }
    bf->Assemble();
    return bf;
}

std::unique_ptr<mfem::ParLinearForm> MechanicalAssembler::AssembleRHS(
    mfem::ParFiniteElementSpace &space, const coeff::MechanicalCoeffs &coeffs,
    BCbuilder::MechanicalBCData &bc)
{
    auto *body_force = coeffs.GetBodyForce();
    auto *lambda = coeffs.GetLambda();
    auto *mu = coeffs.GetMu();
    auto *alpha = coeffs.GetAlpha();
    auto *temperature = coeffs.GetTemperature();
    double reference_temperature = coeffs.GetReferenceTemperature();

    if (!body_force || !lambda || !mu || !alpha || !temperature)
    {
        throw std::runtime_error("MechanicalAssembler::AssembleRHS requires body force, lambda, "
                                 "mu, alpha, and temperature coefficients");
    }

    auto lf = std::make_unique<mfem::ParLinearForm>(&space);
    lf->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(*body_force));
    for (size_t i = 0; i < bc.pressure_markers.size(); ++i)
    {
        lf->AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(bc.pressure_coeffs[i]),
                                  bc.pressure_markers[i]);
    }
    lf->AddDomainIntegrator(new integrator::ThermalStrainLFIntegrator(
        *lambda, *mu, *alpha, *temperature, reference_temperature));
    lf->Assemble();
    return lf;
}

}  // namespace fem::assembler