#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{

struct AssembledSystem
{
    std::unique_ptr<mfem::BilinearForm> bilinear;
    std::unique_ptr<mfem::LinearForm> linear;
    mfem::OperatorPtr A;
    mfem::Vector X;
    mfem::Vector B;
    mfem::Array<int> essential_tdofs;
};

/// Input data for Poisson-type assembly (-div(k grad u) = f).
/// Boundary conditions use Config types directly; the assembler
/// handles marker construction and coefficient creation internally.
struct PoissonAssemblyInput
{
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &diffusion;
    mfem::Coefficient &source;

    // Boundary conditions — Config types, assembler builds markers internally
    std::vector<frontend::ScalarFieldConfig::DirichletBC> dirichlet_bcs;
    std::vector<frontend::ScalarFieldConfig::RobinBC> robin_bcs;

    // Default value for initial guess (e.g. 0.0 for voltage, 293.15 for temperature)
    double default_value = 0.0;

    // Transient: mass-matrix coefficient (ρCp/dt) on LHS.
    mfem::Coefficient *mass_coeff = nullptr;
    // Transient: mass-RHS coefficient (ρCp/dt * T_old) on RHS.
    mfem::Coefficient *mass_rhs_coeff = nullptr;
};

class MfemPoissonAssembler
{
  public:
    static AssembledSystem Assemble(PoissonAssemblyInput &input, mfem::GridFunction &initial_guess);
};

}  // namespace fem::assembly
