#pragma once

#include "spdlog/spdlog.h"

#include <memory>
#include <string>

namespace fem::log
{
/// Initialize the logger with elapsed-time format.
void Init(const std::string &level);

/// Get the global logger instance.
std::shared_ptr<spdlog::logger> Get();
}  // namespace fem::log
