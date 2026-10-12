include_guard(GLOBAL)

include(FetchContent)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # This MFEM revision uses min/max OpenMP reductions, which MSVC's
    # legacy /openmp backend rejects. The compiler remains cl.exe.
    set(OpenMP_RUNTIME_MSVC llvm CACHE STRING "MSVC OpenMP runtime")
    set(OpenMP_CXX_FLAGS "-openmp:llvm" CACHE STRING "MSVC OpenMP C++ flags")
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")

# Keep lightweight libraries in this project's build tree by default.
option(FEM_USE_SYSTEM_LIGHTWEIGHT_DEPENDENCIES
    "Prefer installed spdlog/json/muparser instead of project-local FetchContent" OFF)
if(FEM_USE_SYSTEM_LIGHTWEIGHT_DEPENDENCIES)
    find_package(spdlog CONFIG QUIET)
endif()
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

if(FEM_USE_SYSTEM_LIGHTWEIGHT_DEPENDENCIES)
    find_package(nlohmann_json CONFIG QUIET)
endif()
if(NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

if(FEM_USE_SYSTEM_LIGHTWEIGHT_DEPENDENCIES)
    find_package(muparser CONFIG QUIET)
endif()
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
