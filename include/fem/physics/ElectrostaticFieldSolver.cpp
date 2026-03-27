#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

namespace fem::physics
{
class ElectrostaticFieldSolver
{
  public:
    ElectrostaticFieldSolver(const frontend::ProjectConfig &config);

    // 设置与热场耦合的温度场,并打开耦合标志.电场的耦合指介电常数随温度变化
    void SetCoupledWithThermal(mfem::GridFunction &temperature_field);

    void Solve();

  private:
    std::shared_ptr<frontend::ProjectConfig> config_;
    bool is_coupled_with_thermal_ = false;
    mfem::GridFunction &temperature_field_ = nullptr;
};
