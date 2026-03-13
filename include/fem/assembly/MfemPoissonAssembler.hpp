#pragma once

#include "fem/assembly/IAssembler.hpp"

namespace fem::assembly
{
class MfemPoissonAssembler final : public IAssembler
{
  public:
    AssembledSystem Assemble(PoissonAssemblyInput &input,
                             mfem::GridFunction &initial_guess) override;
};
}  // namespace fem::assembly
