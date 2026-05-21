#include "assembler/MechanicalAssembler.hpp"

#include "assembler/LocalAssemblyExecutor.hpp"
#include "coeff/MechanicalCoeffs.hpp"
#include "integrator/MechanicalIntegrators.hpp"

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

    std::vector<BilinearIntegratorFactory> domain_bf;
    domain_bf.push_back([lambda, mu]()
                        { return std::make_unique<mfem::ElasticityIntegrator>(*lambda, *mu); });

    std::vector<BoundaryBilinearIntegratorSpec> boundary_bf;
    for (size_t i = 0; i < bc.nd_markers.size(); ++i)
    {
        auto *coeff = &bc.penalty_coeffs[i];
        auto *marker = &bc.nd_markers[i];
        boundary_bf.push_back(
            {marker, [coeff]() {
                 return std::make_unique<integrator::CustomNormalDisplacementPenaltyIntegrator>(
                     *coeff);
             }});
    }
    auto bf = AssembleBilinearFormOpenMP(space, domain_bf, boundary_bf);
    bf->Finalize();
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

    std::vector<LinearIntegratorFactory> domain_lf;
    domain_lf.emplace_back([body_force]()
                        { return std::make_unique<mfem::VectorDomainLFIntegrator>(*body_force); });
    domain_lf.emplace_back([lambda, mu, alpha, temperature, reference_temperature]()
        {
            return std::make_unique<integrator::CustomThermalStrainLFIntegrator>(
                *lambda, *mu, *alpha, *temperature, reference_temperature);
        });

    std::vector<BoundaryLinearIntegratorSpec> boundary_lf;
    for (size_t i = 0; i < bc.pressure_markers.size(); ++i)
    {
        auto *coeff = &bc.pressure_coeffs[i];
        auto *marker = &bc.pressure_markers[i];
        boundary_lf.push_back(
            {marker,
             [coeff]() { return std::make_unique<mfem::VectorBoundaryFluxLFIntegrator>(*coeff); }});
    }
    return AssembleLinearFormOpenMP(space, domain_lf, boundary_lf);
}

}  // namespace fem::assembler
