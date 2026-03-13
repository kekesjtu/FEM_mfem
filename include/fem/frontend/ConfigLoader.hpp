#pragma once

#include "fem/frontend/Config.hpp"

#include <string>

namespace fem::frontend
{
class ConfigLoader
{
  public:
    static ProjectConfig LoadFromFile(const std::string &path);
};
}  // namespace fem::frontend
