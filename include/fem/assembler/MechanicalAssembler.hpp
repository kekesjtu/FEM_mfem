#pragma once

#include "fem/BCbuilder/MechanicalBCBuilder.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::coeff
{
class MechanicalCoeffs;
}

namespace fem::assembler
{

class MechanicalAssembler
{
  public:
    static std::unique_ptr<mfem::ParBilinearForm> AssembleStiffness(
        mfem::ParFiniteElementSpace &space, const coeff::MechanicalCoeffs &coeffs,
        BCbuilder::MechanicalBCData &bc);

    static std::unique_ptr<mfem::ParLinearForm> AssembleRHS(mfem::ParFiniteElementSpace &space,
                                                            const coeff::MechanicalCoeffs &coeffs,
                                                            BCbuilder::MechanicalBCData &bc);
};

}  // namespace fem::assembler