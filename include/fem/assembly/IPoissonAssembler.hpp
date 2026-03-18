#pragma once

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{
struct RobinBoundaryCondition
{
    mfem::ConstantCoefficient l;
    mfem::ConstantCoefficient q;
    mfem::Array<int> marker;

    RobinBoundaryCondition(double l_value, double q_value, const mfem::Array<int> &bdr_marker)
        : l(l_value), q(q_value), marker(bdr_marker)
    {
    }
};

struct DirichletBoundaryCondition
{
    mfem::ConstantCoefficient value;
    mfem::Array<int> marker;

    DirichletBoundaryCondition(double value_scalar, const mfem::Array<int> &bdr_marker)
        : value(value_scalar), marker(bdr_marker)
    {
    }
};

struct PoissonAssemblyInput
{
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &diffusion;
    mfem::Coefficient &source;
    mfem::Array<int> essential_bdr_marker;
    std::vector<DirichletBoundaryCondition> dirichlet_conditions;
    std::vector<RobinBoundaryCondition> robin_conditions;
    double dirichlet_value = 0.0;
};

struct AssembledSystem
{
    std::unique_ptr<mfem::BilinearForm> bilinear;
    std::unique_ptr<mfem::LinearForm> linear;
    mfem::OperatorPtr A;
    mfem::Vector X;
    mfem::Vector B;
    mfem::Array<int> essential_tdofs;
};

class IPoissonAssembler
{
  public:
    virtual ~IPoissonAssembler() = default;
    virtual AssembledSystem Assemble(PoissonAssemblyInput &input,
                                     mfem::GridFunction &initial_guess) = 0;
};
}  // namespace fem::assembly
