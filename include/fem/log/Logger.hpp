#pragma once

#include "spdlog/spdlog.h"

#include <memory>
#include <string>

namespace fem::log
{
void Init(const std::string &level);
std::shared_ptr<spdlog::logger> Get();
}  // namespace fem::log
