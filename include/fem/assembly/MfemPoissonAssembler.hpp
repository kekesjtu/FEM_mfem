#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

namespace fem::assembly
{
class MfemPoissonAssembler
{
    struct AssembledSystem
    {
        std::unique_ptr<mfem::BilinearForm> bilinear;
        std::unique_ptr<mfem::LinearForm> linear;
        mfem::OperatorPtr A;
        mfem::Vector X;
        mfem::Vector B;
        mfem::Array<int> essential_tdofs;
    };

  public:
    MfemPoissonAssembler(frontend::ProjectConfig &config)
        : config_(std::make_shared<frontend::ProjectConfig>(config)){};
    AssembledSystem Assemble();

  private:
    std::shared_ptr<frontend::ProjectConfig> config_;
};
}  // namespace fem::assembly
