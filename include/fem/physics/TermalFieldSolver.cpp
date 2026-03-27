#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

namespace fem::physics
{
class ThermalFieldSolver
{
  public:
    ThermalFieldSolver(const frontend::ProjectConfig &config);

    // 设置与电场耦合的温度场,并打开耦合标志.热场的耦合指源项用电场的焦耳热代替
    void SetCoupledWithElectrostatic(mfem::GridFunction &voltage_field);

    void Solve();

  private:
    void CalculateJouleHeating();  // 由电场计算焦耳热源项,并更新线性形式

  private:
    std::shared_ptr<frontend::ProjectConfig> config_;
    bool is_coupled_with_electrostatic_ = false;
    mfem::GridFunction &voltage_field_ = nullptr;
    mfem::LinearForm &joule_heating_ = nullptr; // 用于更新焦耳热源项
};
