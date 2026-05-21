#pragma once

#include "mfem.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace fem::assembler
{

using BilinearIntegratorFactory =
    std::function<std::unique_ptr<mfem::BilinearFormIntegrator>()>;
using LinearIntegratorFactory = std::function<std::unique_ptr<mfem::LinearFormIntegrator>()>;

struct BoundaryBilinearIntegratorSpec
{
    const mfem::Array<int> *marker = nullptr;
    BilinearIntegratorFactory factory;
};

struct BoundaryLinearIntegratorSpec
{
    const mfem::Array<int> *marker = nullptr;
    LinearIntegratorFactory factory;
};

std::unique_ptr<mfem::ParBilinearForm> AssembleBilinearFormOpenMP(
    mfem::ParFiniteElementSpace &space,
    const std::vector<BilinearIntegratorFactory> &domain_integrators,
    const std::vector<BoundaryBilinearIntegratorSpec> &boundary_integrators);

std::unique_ptr<mfem::ParLinearForm> AssembleLinearFormOpenMP(
    mfem::ParFiniteElementSpace &space, const std::vector<LinearIntegratorFactory> &domain_integrators,
    const std::vector<BoundaryLinearIntegratorSpec> &boundary_integrators);

}  // namespace fem::assembler
