#include "fem/assembly/MechanicalAssembler.hpp"

#include "fem/integrator/MechanicalIntegrators.hpp"

namespace fem::assembly
{

std::unique_ptr<mfem::ParBilinearForm> MechanicalAssembler::AssembleStiffness(
    mfem::ParFiniteElementSpace &space, mfem::Coefficient &lambda, mfem::Coefficient &mu,
    std::vector<mfem::ConstantCoefficient> &penalty_coeffs,
    std::vector<mfem::Array<int>> &nd_markers)
{
    auto bf = std::make_unique<mfem::ParBilinearForm>(&space);
    bf->AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
    for (size_t i = 0; i < nd_markers.size(); ++i)
    {
        bf->AddBoundaryIntegrator(
            new integrator::NormalDisplacementPenaltyIntegrator(penalty_coeffs[i]), nd_markers[i]);
    }
    bf->Assemble();
    return bf;
}

std::unique_ptr<mfem::ParLinearForm> MechanicalAssembler::AssembleRHS(
    mfem::ParFiniteElementSpace &space, mfem::VectorCoefficient &body_force,
    std::vector<mfem::ConstantCoefficient> &pressure_coeffs,
    std::vector<mfem::Array<int>> &pressure_markers, mfem::Coefficient *lambda,
    mfem::Coefficient *mu, mfem::Coefficient *alpha, mfem::Coefficient *temperature,
    double reference_temperature)
{
    auto lf = std::make_unique<mfem::ParLinearForm>(&space);
    lf->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(body_force));
    for (size_t i = 0; i < pressure_markers.size(); ++i)
    {
        lf->AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(pressure_coeffs[i]),
                                  pressure_markers[i]);
    }
    if (alpha && temperature && lambda && mu)
    {
        lf->AddDomainIntegrator(new integrator::ThermalStrainLFIntegrator(
            *lambda, *mu, *alpha, *temperature, reference_temperature));
    }
    lf->Assemble();
    return lf;
}

}  // namespace fem::assembly
