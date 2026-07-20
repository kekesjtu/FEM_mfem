include_guard(GLOBAL)

include(FetchContent)

find_package(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")

# Prefer installed Config packages. Fall back to pinned upstream versions when
# a package is not installed on the host.
find_package(spdlog CONFIG QUIET)
if(NOT TARGET spdlog::spdlog)
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "Disable spdlog examples" FORCE)
    set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "Disable spdlog tests" FORCE)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.15.3
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(spdlog)
endif()

find_package(nlohmann_json CONFIG QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

find_package(muparser CONFIG QUIET)
if(NOT TARGET muparser::muparser)
    set(ENABLE_SAMPLES OFF CACHE BOOL "Disable muparser samples" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests" FORCE)
    FetchContent_Declare(
        muparser
        GIT_REPOSITORY https://github.com/beltoforion/muparser.git
        GIT_TAG v2.3.5
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(muparser)
endif()
