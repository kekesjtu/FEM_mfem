#pragma once

#include "fem/assembly/IMechanicalAssembler.hpp"

namespace fem::assembly
{
class MfemLinearElasticityAssembler : public IMechanicalAssembler
{
  public:
    MechanicalSystem Assemble(LinearElasticityInput &input,
                              mfem::GridFunction &initial_guess) override;
};
}  // namespace fem::assembly
