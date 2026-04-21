#pragma once

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace fem::assembly
{

class MechanicalAssembler
{
  public:
    static std::unique_ptr<mfem::ParBilinearForm> AssembleStiffness(
        mfem::ParFiniteElementSpace &space, mfem::Coefficient &lambda, mfem::Coefficient &mu,
        std::vector<mfem::ConstantCoefficient> &penalty_coeffs,
        std::vector<mfem::Array<int>> &nd_markers);

    static std::unique_ptr<mfem::ParLinearForm> AssembleRHS(
        mfem::ParFiniteElementSpace &space, mfem::VectorCoefficient &body_force,
        std::vector<mfem::ConstantCoefficient> &pressure_coeffs,
        std::vector<mfem::Array<int>> &pressure_markers, mfem::Coefficient *lambda,
        mfem::Coefficient *mu, mfem::Coefficient *alpha, mfem::Coefficient *temperature,
        double reference_temperature);
};

}  // namespace fem::assembly
