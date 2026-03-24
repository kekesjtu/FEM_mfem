#pragma once

#include "fem/assembly/IPoissonAssembler.hpp"

namespace fem::assembly
{
class MfemParallelPoissonAssembler final : public IPoissonAssembler
{
  public:
    AssembledSystem Assemble(PoissonAssemblyInput &input,
                             mfem::GridFunction &initial_guess) override;
};
}  // namespace fem::assembly
