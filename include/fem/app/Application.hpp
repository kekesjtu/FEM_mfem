#pragma once

#include <string>

namespace fem::app
{
class Application
{
  public:
    int Run(const std::string &config_path);
};
}  // namespace fem::app
