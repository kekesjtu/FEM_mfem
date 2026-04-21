#pragma once

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{

struct AssembledSystem
{
    std::unique_ptr<mfem::ParBilinearForm> bilinear;
    std::unique_ptr<mfem::ParLinearForm> linear;
    mfem::OperatorPtr A;
    mfem::Vector X;
    mfem::Vector B;
    mfem::Array<int> essential_tdofs;
};

struct ElectrostaticBCData
{
    mfem::Array<int> essential_bdr;
    std::vector<mfem::ConstantCoefficient> dirichlet_coeffs;
    std::vector<mfem::Array<int>> dirichlet_markers;
    std::vector<mfem::ConstantCoefficient> robin_l_coeffs;
    std::vector<mfem::ConstantCoefficient> robin_q_coeffs;
    std::vector<mfem::Array<int>> robin_markers;
};

class ElectrostaticAssembler
{
  public:
    static AssembledSystem Assemble(mfem::ParFiniteElementSpace &space, mfem::Coefficient &sigma,
                                    mfem::Coefficient &source, ElectrostaticBCData &bc,
                                    mfem::GridFunction &voltage);
};

}  // namespace fem::assembly
