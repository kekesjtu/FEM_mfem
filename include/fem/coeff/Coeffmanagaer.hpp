// 这个文件用于将config中的material
#pragma once

#include "fem/frontend/Config.hpp"
#include "mfem.hpp"

namespace fem::coeff
{

class CoefficientManager
{
  public:
    // 根据配置文件中的材料信息创建相应的系数对象
    // 分域常数系数 -用于一般的场求解
    // 分域表达式/常数系数 -用于电热耦合中,电导率随温度变化

};

}  // namespace fem::coeff