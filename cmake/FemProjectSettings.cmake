include_guard(GLOBAL)

# Generate compile_commands.json for editors and static-analysis tools.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(FEM_DEPENDENCY_ROOT "$ENV{FEM_DEPENDENCY_ROOT}" CACHE PATH
    "Native dependency install prefix, or parent of individual package prefixes")
if(FEM_DEPENDENCY_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH "${FEM_DEPENDENCY_ROOT}")
    foreach(_dependency IN ITEMS hypre metis mumps)
        if(IS_DIRECTORY "${FEM_DEPENDENCY_ROOT}/${_dependency}")
            list(APPEND CMAKE_PREFIX_PATH "${FEM_DEPENDENCY_ROOT}/${_dependency}")
        endif()
    endforeach()
endif()

# Keep project warnings target-scoped so they are not imposed on MFEM or
# dependencies brought in through FetchContent.
add_library(fem_project_options INTERFACE)
if(MSVC)
    target_compile_options(fem_project_options INTERFACE /W4 /utf-8 /permissive-)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fem_project_options INTERFACE -Wall -fdiagnostics-color=auto)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(fem_project_options INTERFACE -Wall -fcolor-diagnostics)
endif()
if(WIN32)
    # Expression.cpp uses M_PI; Windows headers also define min/max macros.
    target_compile_definitions(fem_project_options INTERFACE _USE_MATH_DEFINES NOMINMAX)
endif()

function(fem_print_configuration)
    message(STATUS "Platform/compiler: ${CMAKE_SYSTEM_NAME} / ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "Dependency root: ${FEM_DEPENDENCY_ROOT}")
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
    message(STATUS "Source directory: ${PROJECT_SOURCE_DIR}")
    message(STATUS "Binary directory: ${PROJECT_BINARY_DIR}")
    message(STATUS "MUMPS backend: ${FEM_ENABLE_MUMPS}")
    message(STATUS "Cluster PARDISO backend: ${FEM_ENABLE_CPARDISO}")
endfunction()
