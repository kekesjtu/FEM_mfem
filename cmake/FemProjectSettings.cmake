include_guard(GLOBAL)

# Generate compile_commands.json for editors and static-analysis tools.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Keep project warnings target-scoped so they are not imposed on MFEM or
# dependencies brought in through FetchContent.
add_library(fem_project_options INTERFACE)
target_compile_options(fem_project_options INTERFACE
    -Wall
    -fdiagnostics-color=auto
)

function(fem_print_configuration)
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
    message(STATUS "Source directory: ${PROJECT_SOURCE_DIR}")
    message(STATUS "Binary directory: ${PROJECT_BINARY_DIR}")
    message(STATUS "MUMPS backend: ${FEM_ENABLE_MUMPS}")
    message(STATUS "Cluster PARDISO backend: ${FEM_ENABLE_CPARDISO}")
endfunction()
