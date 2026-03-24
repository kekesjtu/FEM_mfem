#pragma once

#include "fem/frontend/Expression.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::mapping
{
class CoefficientFactory
{
  public:
    static std::unique_ptr<mfem::Coefficient> CreateScalar(
        const fem::frontend::ScalarExpression &expr);
};
}  // namespace fem::mapping
