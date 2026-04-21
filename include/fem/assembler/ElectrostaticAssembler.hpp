#pragma once

#include "fem/BCbuilder/ScalarBCBuilder.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::coeff
{
class ElectrostaticCoeffs;
}

namespace fem::assembler
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

class ElectrostaticAssembler
{
  public:
    static AssembledSystem Assemble(mfem::ParFiniteElementSpace &space,
                                    coeff::ElectrostaticCoeffs &coeffs, BCbuilder::ScalarBCData &bc,
                                    mfem::GridFunction &voltage);
};

}  // namespace fem::assembler