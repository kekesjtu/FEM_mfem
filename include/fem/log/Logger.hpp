#pragma once

#include <memory>
#include <string>

namespace spdlog
{
class logger;
}

namespace fem::log
{
void Init(const std::string &level);
std::shared_ptr<spdlog::logger> Get();
}  // namespace fem::log
