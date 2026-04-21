#pragma once

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
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
                                            mfem::Coefficient &k, mfem::Coefficient &rho_cp,
                                            std::vector<mfem::ConstantCoefficient> &robin_l_coeffs,
                                            std::vector<mfem::Array<int>> &robin_markers);

    static mfem::Vector AssembleSourceRHS(mfem::ParFiniteElementSpace &fespace,
                                          mfem::Coefficient &source,
                                          std::vector<mfem::ConstantCoefficient> &robin_q_coeffs,
                                          std::vector<mfem::Array<int>> &robin_markers);
};

}  // namespace fem::assembly
