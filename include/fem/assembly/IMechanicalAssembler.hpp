#pragma once

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{
struct MechanicalTractionBoundary
{
    mfem::VectorConstantCoefficient traction;
    mfem::Array<int> marker;

    MechanicalTractionBoundary(const mfem::Vector &traction_value,
                               const mfem::Array<int> &bdr_marker)
        : traction(traction_value), marker(bdr_marker)
    {
    }
};

struct MechanicalPressureBoundary
{
    mfem::ConstantCoefficient pressure;
    mfem::Array<int> marker;

    MechanicalPressureBoundary(double pressure_value, const mfem::Array<int> &bdr_marker)
        : pressure(pressure_value), marker(bdr_marker)
    {
    }
};

struct MechanicalNormalDisplacementBoundary
{
    mfem::ConstantCoefficient penalty;
    mfem::Array<int> marker;

    MechanicalNormalDisplacementBoundary(double penalty_value, const mfem::Array<int> &bdr_marker)
        : penalty(penalty_value), marker(bdr_marker)
    {
    }
};

struct LinearElasticityInput
{
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &lambda;
    mfem::Coefficient &mu;
    mfem::VectorCoefficient &body_force;
    mfem::Array<int> essential_bdr_marker;
    mfem::VectorConstantCoefficient &dirichlet_displacement;
    std::vector<MechanicalTractionBoundary> traction_boundaries;
    std::vector<MechanicalPressureBoundary> pressure_boundaries;
    std::vector<MechanicalNormalDisplacementBoundary> normal_displacement_boundaries;
    mfem::Coefficient *thermal_expansion_secant = nullptr;
    mfem::Coefficient *temperature = nullptr;
    double reference_temperature = 293.15;
};

struct MechanicalSystem
{
    std::unique_ptr<mfem::BilinearForm> bilinear;
    std::unique_ptr<mfem::LinearForm> linear;
    mfem::OperatorPtr A;
    mfem::Vector X;
    mfem::Vector B;
    mfem::Array<int> essential_tdofs;
};

class IMechanicalAssembler
{
  public:
    virtual ~IMechanicalAssembler() = default;

    virtual MechanicalSystem Assemble(LinearElasticityInput &input,
                                      mfem::GridFunction &initial_guess) = 0;
};
}  // namespace fem::assembly
