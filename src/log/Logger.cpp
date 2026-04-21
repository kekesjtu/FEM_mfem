#include "fem/log/Logger.hpp"

#include "mfem.hpp"
#include "spdlog/pattern_formatter.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <algorithm>
#include <chrono>

namespace fem::log
{
namespace
{

/// Global program start time for elapsed-time logging.
const auto kProgramStart = std::chrono::steady_clock::now();

std::string FormatElapsedTime(std::chrono::steady_clock::duration elapsed)
{
    const double secs_f = std::chrono::duration<double>(elapsed).count();
    return fmt::format("{:.3f}s", secs_f);
}

/// Custom spdlog flag formatter: outputs elapsed seconds since program start.
/// Usage in pattern: %* -> "4.877s", "1:04.877", "2:01:05.012"
class ElapsedTimeFlag : public spdlog::custom_flag_formatter
{
  public:
    void format(const spdlog::details::log_msg & /*msg*/, const std::tm & /*tm*/,
                spdlog::memory_buf_t &dest) override
    {
        auto elapsed = std::chrono::steady_clock::now() - kProgramStart;
        auto txt = FormatElapsedTime(elapsed);
        dest.append(txt.data(), txt.data() + txt.size());
    }

    std::unique_ptr<custom_flag_formatter> clone() const override
    {
        return std::make_unique<ElapsedTimeFlag>();
    }
};

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
        return spdlog::level::trace;
    if (level == "debug")
        return spdlog::level::debug;
    if (level == "warn")
        return spdlog::level::warn;
    if (level == "error")
        return spdlog::level::err;
    if (level == "critical")
        return spdlog::level::critical;
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

    auto effective_level = ParseLevel(level);

    // Suppress logging on non-root ranks
    if (mfem::Mpi::IsInitialized() && mfem::Mpi::WorldRank() != 0)
    {
        effective_level = spdlog::level::off;
    }

    logger->set_level(effective_level);

    // Elapsed-time pattern: [4.877s] [info] message
    auto formatter = std::make_unique<spdlog::pattern_formatter>();
    formatter->add_flag<ElapsedTimeFlag>('*');
    formatter->set_pattern("[%*] [%^%l%$] %v");
    logger->set_formatter(std::move(formatter));

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
