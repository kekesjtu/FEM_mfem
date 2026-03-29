#pragma once

#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{

/// Input data for linear elasticity assembly.
/// Boundary conditions use Config types directly; the assembler
/// handles marker construction and coefficient creation internally.
struct LinearElasticityInput
{
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &lambda;
    mfem::Coefficient &mu;
    mfem::VectorCoefficient &body_force;

    // Optional thermal coupling
    mfem::Coefficient *thermal_expansion_secant = nullptr;
    mfem::Coefficient *temperature = nullptr;
    double reference_temperature = 293.15;

    // Boundary conditions — Config types, assembler builds markers internally
    std::vector<frontend::MechanicalFieldConfig::DisplacementBC> displacement_bcs;
    std::vector<frontend::MechanicalFieldConfig::NormalDisplacementBC> normal_displacement_bcs;
    std::vector<frontend::MechanicalFieldConfig::PressureBC> pressure_bcs;
};

/// Assembled system for mechanical problems (reuses AssembledSystem from Poisson)
using MechanicalSystem = AssembledSystem;

class MfemLinearElasticityAssembler
{
  public:
    static MechanicalSystem Assemble(LinearElasticityInput &input,
                                     mfem::GridFunction &initial_guess);
};

}  // namespace fem::assembly
