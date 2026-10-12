cmake_minimum_required(VERSION 3.23)

# Use the DLL CRT consistently with native Windows dependencies (/MD, /MDd).
# This must be initialized before project() enables the compilers.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
    CACHE STRING "MSVC runtime library")

project(FEM_mfem
    VERSION 1.0.0
    DESCRIPTION "Finite Element Method Solver based on MFEM"
    LANGUAGES C CXX
)

cmake_path(IS_PREFIX PROJECT_SOURCE_DIR "${PROJECT_BINARY_DIR}" NORMALIZE _fem_build_inside_project)
if(NOT _fem_build_inside_project OR PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR
        "Build this project in its own build directory, for example "
        "cmake -S . -B build/windows-msvc-intel-full-release. "
        "The shared dependency directory must not contain application build outputs.")
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FemProjectSettings)
include(FemDependencies)
include(FemMfem)
include(FemTarget)

fem_print_configuration()
