#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/integrator/PoissonIntegrators.hpp"

namespace fem::assembly
{
AssembledSystem MfemPoissonAssembler::Assemble(PoissonAssemblyInput &input,
                                               mfem::GridFunction &initial_guess)
{
    AssembledSystem system;
    const int num_bdr = input.space.GetMesh()->bdr_attributes.Max();

    // --- Build Robin BC markers and coefficients from Config types ---
    std::vector<mfem::ConstantCoefficient> robin_l_coeffs, robin_q_coeffs;
    std::vector<mfem::Array<int>> robin_markers;
    robin_l_coeffs.reserve(input.robin_bcs.size());
    robin_q_coeffs.reserve(input.robin_bcs.size());
    robin_markers.reserve(input.robin_bcs.size());

    for (const auto &rbc : input.robin_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : rbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
            {
                marker[attr - 1] = 1;
            }
        }
        robin_l_coeffs.emplace_back(rbc.l);
        robin_q_coeffs.emplace_back(rbc.q);
        robin_markers.push_back(std::move(marker));
    }

    // --- Build Dirichlet BC markers and coefficients ---
    mfem::Array<int> essential_bdr(num_bdr);
    essential_bdr = 0;
    std::vector<mfem::ConstantCoefficient> dirichlet_coeffs;
    std::vector<mfem::Array<int>> dirichlet_markers;
    dirichlet_coeffs.reserve(input.dirichlet_bcs.size());
    dirichlet_markers.reserve(input.dirichlet_bcs.size());

    for (const auto &dbc : input.dirichlet_bcs)
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
        dirichlet_coeffs.emplace_back(dbc.value);
        dirichlet_markers.push_back(std::move(marker));
    }

#ifdef MFEM_USE_MPI
    auto *par_fespace = dynamic_cast<mfem::ParFiniteElementSpace *>(&input.space);
    if (par_fespace)
    {
        // --- Parallel assembly ---
        auto par_linear = std::make_unique<mfem::ParLinearForm>(par_fespace);
        par_linear->AddDomainIntegrator(new integrator::CustomDomainLFIntegrator(input.source));
        if (input.mass_rhs_coeff)
        {
            par_linear->AddDomainIntegrator(new mfem::DomainLFIntegrator(*input.mass_rhs_coeff));
        }
        for (size_t i = 0; i < robin_markers.size(); ++i)
        {
            par_linear->AddBoundaryIntegrator(
                new integrator::CustomBoundaryLFIntegrator(robin_q_coeffs[i]), robin_markers[i]);
        }
        par_linear->Assemble();

        auto par_bilinear = std::make_unique<mfem::ParBilinearForm>(par_fespace);
        par_bilinear->AddDomainIntegrator(
            new integrator::CustomDiffusionIntegrator(input.diffusion));
        if (input.mass_coeff)
        {
            par_bilinear->AddDomainIntegrator(new mfem::MassIntegrator(*input.mass_coeff));
        }
        for (size_t i = 0; i < robin_markers.size(); ++i)
        {
            par_bilinear->AddBoundaryIntegrator(
                new integrator::CustomBoundaryMassIntegrator(robin_l_coeffs[i]), robin_markers[i]);
        }
        par_bilinear->Assemble();

        par_fespace->GetEssentialTrueDofs(essential_bdr, system.essential_tdofs);
        initial_guess = input.default_value;
        for (size_t i = 0; i < dirichlet_markers.size(); ++i)
        {
            initial_guess.ProjectBdrCoefficient(dirichlet_coeffs[i], dirichlet_markers[i]);
        }

        par_bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *par_linear, system.A,
                                       system.X, system.B);

        system.bilinear = std::move(par_bilinear);
        system.linear = std::move(par_linear);
        return system;
    }
#endif

    // --- Serial assembly ---
    system.linear = std::make_unique<mfem::LinearForm>(&input.space);
    system.linear->AddDomainIntegrator(new integrator::CustomDomainLFIntegrator(input.source));
    if (input.mass_rhs_coeff)
    {
        system.linear->AddDomainIntegrator(new mfem::DomainLFIntegrator(*input.mass_rhs_coeff));
    }
    for (size_t i = 0; i < robin_markers.size(); ++i)
    {
        system.linear->AddBoundaryIntegrator(
            new integrator::CustomBoundaryLFIntegrator(robin_q_coeffs[i]), robin_markers[i]);
    }
    system.linear->Assemble();

    system.bilinear = std::make_unique<mfem::BilinearForm>(&input.space);
    system.bilinear->AddDomainIntegrator(
        new integrator::CustomDiffusionIntegrator(input.diffusion));
    if (input.mass_coeff)
    {
        system.bilinear->AddDomainIntegrator(new mfem::MassIntegrator(*input.mass_coeff));
    }
    for (size_t i = 0; i < robin_markers.size(); ++i)
    {
        system.bilinear->AddBoundaryIntegrator(
            new integrator::CustomBoundaryMassIntegrator(robin_l_coeffs[i]), robin_markers[i]);
    }
    system.bilinear->Assemble();

    input.space.GetEssentialTrueDofs(essential_bdr, system.essential_tdofs);
    initial_guess = input.default_value;
    for (size_t i = 0; i < dirichlet_markers.size(); ++i)
    {
        initial_guess.ProjectBdrCoefficient(dirichlet_coeffs[i], dirichlet_markers[i]);
    }

    system.bilinear->FormLinearSystem(system.essential_tdofs, initial_guess, *system.linear,
                                      system.A, system.X, system.B);

    return system;
}
}  // namespace fem::assembly
