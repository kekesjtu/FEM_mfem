#pragma once

#include "mfem.hpp"

#include <memory>

namespace fem::assembly
{
struct PoissonAssemblyInput
{
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &diffusion;
    mfem::Coefficient &source;
    mfem::Array<int> essential_bdr_marker;
    mfem::Coefficient *robin_l = nullptr;
    mfem::Coefficient *robin_q = nullptr;
    mfem::Array<int> robin_bdr_marker;
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

class IAssembler
{
  public:
    virtual ~IAssembler() = default;
    virtual AssembledSystem Assemble(PoissonAssemblyInput &input,
                                     mfem::GridFunction &initial_guess) = 0;
};
}  // namespace fem::assembly
