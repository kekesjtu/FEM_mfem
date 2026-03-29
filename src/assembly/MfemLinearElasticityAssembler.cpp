#include "fem/assembly/MfemLinearElasticityAssembler.hpp"
#include "fem/integrator/LinearElasticityIntegrators.hpp"

namespace fem::assembly
{

namespace
{

// Build normal displacement penalty markers and coefficients from Config types
void BuildNormalDisplacementData(
    const std::vector<frontend::MechanicalFieldConfig::NormalDisplacementBC> &ndbcs, int num_bdr,
    std::vector<mfem::ConstantCoefficient> &penalty_coeffs,
    std::vector<mfem::Array<int>> &nd_markers)
{
    penalty_coeffs.reserve(ndbcs.size());
    nd_markers.reserve(ndbcs.size());
    for (const auto &nbc : ndbcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : nbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        penalty_coeffs.emplace_back(nbc.penalty);
        nd_markers.push_back(std::move(marker));
    }
}

// Build pressure boundary markers and coefficients
void BuildPressureData(const std::vector<frontend::MechanicalFieldConfig::PressureBC> &pbcs,
                       int num_bdr, std::vector<mfem::ConstantCoefficient> &pressure_coeffs,
                       std::vector<mfem::Array<int>> &pressure_markers)
{
    pressure_coeffs.reserve(pbcs.size());
    pressure_markers.reserve(pbcs.size());
    for (const auto &pbc : pbcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : pbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        pressure_coeffs.emplace_back(pbc.value);
        pressure_markers.push_back(std::move(marker));
    }
}

// Build displacement essential BC markers and coefficients
void BuildDisplacementData(const std::vector<frontend::MechanicalFieldConfig::DisplacementBC> &dbcs,
                           int num_bdr, int dim, mfem::Array<int> &essential_bdr,
                           std::vector<mfem::VectorConstantCoefficient> &disp_coeffs,
                           std::vector<mfem::Array<int>> &disp_markers)
{
    disp_coeffs.reserve(dbcs.size());
    disp_markers.reserve(dbcs.size());
    for (const auto &dbc : dbcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : dbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
            {
                marker[attr - 1] = 1;
                essential_bdr[attr - 1] = 1;
            }
        }
        mfem::Vector val(dim);
        for (int d = 0; d < dim && d < static_cast<int>(dbc.value.size()); ++d)
            val(d) = dbc.value[d];
        disp_coeffs.emplace_back(val);
        disp_markers.push_back(std::move(marker));
    }
}

// Add LHS integrators to a bilinear form (works for both serial and parallel)
void AddLHSIntegrators(mfem::BilinearForm &bf, LinearElasticityInput &input,
                       std::vector<mfem::ConstantCoefficient> &penalty_coeffs,
                       std::vector<mfem::Array<int>> &nd_markers)
{
    bf.AddDomainIntegrator(new mfem::ElasticityIntegrator(input.lambda, input.mu));
    for (size_t i = 0; i < nd_markers.size(); ++i)
    {
        bf.AddBoundaryIntegrator(
            new integrator::NormalDisplacementPenaltyIntegrator(penalty_coeffs[i]), nd_markers[i]);
    }
    bf.Assemble();
}

// Add RHS integrators to a linear form (works for both serial and parallel)
void AddRHSIntegrators(mfem::LinearForm &lf, LinearElasticityInput &input,
                       std::vector<mfem::ConstantCoefficient> &pressure_coeffs,
                       std::vector<mfem::Array<int>> &pressure_markers)
{
    lf.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(input.body_force));
    for (size_t i = 0; i < pressure_markers.size(); ++i)
    {
        lf.AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(pressure_coeffs[i]),
                                 pressure_markers[i]);
    }
    if (input.thermal_expansion_secant && input.temperature)
    {
        lf.AddDomainIntegrator(new integrator::ThermalStrainLFIntegrator(
            input.lambda, input.mu, *input.thermal_expansion_secant, *input.temperature,
            input.reference_temperature));
    }
    lf.Assemble();
}

}  // namespace

MechanicalSystem MfemLinearElasticityAssembler::Assemble(LinearElasticityInput &input,
                                                         mfem::GridFunction &initial_guess)
{
    AssembledSystem system;
    const int num_bdr = input.space.GetMesh()->bdr_attributes.Max();
    const int dim = input.space.GetMesh()->Dimension();

    // --- Build BC data ---
    std::vector<mfem::ConstantCoefficient> penalty_coeffs;
    std::vector<mfem::Array<int>> nd_markers;
    BuildNormalDisplacementData(input.normal_displacement_bcs, num_bdr, penalty_coeffs, nd_markers);

    std::vector<mfem::ConstantCoefficient> pressure_coeffs;
    std::vector<mfem::Array<int>> pressure_markers;
    BuildPressureData(input.pressure_bcs, num_bdr, pressure_coeffs, pressure_markers);

    mfem::Array<int> essential_bdr(num_bdr);
    essential_bdr = 0;
    std::vector<mfem::VectorConstantCoefficient> disp_coeffs;
    std::vector<mfem::Array<int>> disp_markers;
    BuildDisplacementData(input.displacement_bcs, num_bdr, dim, essential_bdr, disp_coeffs,
                          disp_markers);

#ifdef MFEM_USE_MPI
    auto *par_fespace = dynamic_cast<mfem::ParFiniteElementSpace *>(&input.space);
    if (par_fespace)
    {
        auto par_bilinear = std::make_unique<mfem::ParBilinearForm>(par_fespace);
        AddLHSIntegrators(*par_bilinear, input, penalty_coeffs, nd_markers);

        auto par_linear = std::make_unique<mfem::ParLinearForm>(par_fespace);
        AddRHSIntegrators(*par_linear, input, pressure_coeffs, pressure_markers);

        initial_guess = 0.0;
        for (size_t i = 0; i < disp_markers.size(); ++i)
            initial_guess.ProjectBdrCoefficient(disp_coeffs[i], disp_markers[i]);

        par_fespace->GetEssentialTrueDofs(essential_bdr, system.essential_tdofs);
        par_bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *par_linear, system.A,
                                       system.X, system.B);

        system.bilinear = std::move(par_bilinear);
        system.linear = std::move(par_linear);
        return system;
    }
#endif

    // --- Serial path ---
    system.bilinear = std::make_unique<mfem::BilinearForm>(&input.space);
    AddLHSIntegrators(*system.bilinear, input, penalty_coeffs, nd_markers);

    system.linear = std::make_unique<mfem::LinearForm>(&input.space);
    AddRHSIntegrators(*system.linear, input, pressure_coeffs, pressure_markers);

    initial_guess = 0.0;
    for (size_t i = 0; i < disp_markers.size(); ++i)
        initial_guess.ProjectBdrCoefficient(disp_coeffs[i], disp_markers[i]);

    input.space.GetEssentialTrueDofs(essential_bdr, system.essential_tdofs);
    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
}
}  // namespace fem::assembly
