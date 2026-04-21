#pragma once

#include "fem/frontend/Config.hpp"

namespace fem::steady
{

class SteadySolver
{
  public:
    explicit SteadySolver(frontend::ProjectConfig& config);

    void SolveSteady();

  private:
    frontend::ProjectConfig& config_;
};

}  // namespace fem::steady
