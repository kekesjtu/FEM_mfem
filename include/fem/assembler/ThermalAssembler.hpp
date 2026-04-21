#pragma once

#include "fem/BCbuilder/ScalarBCBuilder.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::coeff
{
class ThermalCoeffs;
}

namespace fem::assembler
{

struct ThermalMatrices
{
    std::unique_ptr<mfem::ParBilinearForm> K_bf;
    std::unique_ptr<mfem::ParBilinearForm> C_bf;
    std::unique_ptr<mfem::HypreParMatrix> K;
    std::unique_ptr<mfem::HypreParMatrix> C;
};

class ThermalAssembler
{
  public:
    static ThermalMatrices AssembleMatrices(mfem::ParFiniteElementSpace &fespace,
                                            const coeff::ThermalCoeffs &coeffs,
                                            BCbuilder::ScalarBCData &bc);

    static mfem::Vector AssembleSourceRHS(mfem::ParFiniteElementSpace &fespace,
                                          const coeff::ThermalCoeffs &coeffs,
                                          BCbuilder::ScalarBCData &bc);
};

}  // namespace fem::assembler
