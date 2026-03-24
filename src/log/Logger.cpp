#include "fem/log/Logger.hpp"

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <algorithm>

namespace fem::log
{
namespace
{
std::shared_ptr<spdlog::logger> &LoggerRef()
{
    static std::shared_ptr<spdlog::logger> logger;
    return logger;
}

spdlog::level::level_enum ParseLevel(std::string level)
{
    std::transform(level.begin(), level.end(), level.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (level == "trace")
    {
        return spdlog::level::trace;
    }
    if (level == "debug")
    {
        return spdlog::level::debug;
    }
    if (level == "warn")
    {
        return spdlog::level::warn;
    }
    if (level == "error")
    {
        return spdlog::level::err;
    }
    if (level == "critical")
    {
        return spdlog::level::critical;
    }
    return spdlog::level::info;
}
}  // namespace

void Init(const std::string &level)
{
    auto &logger = LoggerRef();
    if (!logger)
    {
        logger = spdlog::stdout_color_mt("fem");
    }

    logger->set_level(ParseLevel(level));
    logger->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");
    spdlog::set_default_logger(logger);
}

std::shared_ptr<spdlog::logger> Get()
{
    auto &logger = LoggerRef();
    if (!logger)
    {
        Init("info");
    }
    return logger;
}
}  // namespace fem::log
